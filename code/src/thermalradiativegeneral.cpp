#include "thermalradiativegeneral.h"
#include "quadraturegrid.h"

ThermalRadiativeGeneral::ThermalRadiativeGeneral( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    _settings->SetSource( true );
    _suOlson         = false;
    _constitutiveLaw = 2;    // 1 is Su Olson, 2 is constant

    // physical constants
    _c             = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a             = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef          = 1.0;                    // reference temperature
    _sigma         = 1.0;                    // opacity
    _alpha         = 4.0 * _a;               // closure relation
    _cV            = 0.718 * 1e7;            // heat capacity. Air has cV = 0.718 in [kJ/(kg K)] = [1000 m^2 / s^2 / K]
    double sigmaSB = 5.6704 * 1e-5;          // Stefan Boltzmann constant in [erg/cm^2/s/K^4]
    _a             = 4.0 * sigmaSB / _c;

    if( !_suOlson ) _TRef = pow( _cV / _a, 1.0 / 4.0 );    // ensure eTilde = O(1)

    _epsilon = 4.0 * _a / _alpha;

    // compute xi Quadrature points
    Vector xiEta( _settings->GetNDimXi() );
    _variances = _settings->GetSigma();

    // get quadrature grid
    auto grid = QuadratureGrid::Create( _settings, _settings->GetNQuadPoints() );
    _xiQuad   = grid->GetNodes();

    // compute Roe flux components
    Matrix P( 2, 2, 1.0 );
    P( 0, 0 ) = -sqrt( 3 ) / _c;
    P( 0, 1 ) = sqrt( 3 ) / _c;
    Matrix PInv( 2, 2, 0.5 );
    PInv( 0, 0 ) = -_c / sqrt( 3 ) / 2.0;
    PInv( 1, 0 ) = _c / sqrt( 3 ) / 2.0;
    Matrix LambdaAbs( 2, 2, 0.0 );
    LambdaAbs( 0, 0 ) = fabs( -1.0 / sqrt( 3 ) / _epsilon );
    LambdaAbs( 1, 1 ) = fabs( 1.0 / sqrt( 3 ) / _epsilon );
    Matrix AbsAPart   = PInv * LambdaAbs * P;
    _AbsA             = Matrix( _nStates, _nStates, 0.0 );
    _AbsA( 0, 0 )     = AbsAPart( 0, 0 );
    _AbsA( 0, 1 )     = AbsAPart( 0, 1 );
    _AbsA( 1, 0 )     = AbsAPart( 1, 0 );
    _AbsA( 1, 1 )     = AbsAPart( 1, 1 );

    // scale refinement threshholds
    _settings->SetRefinementThreshold( 1e-15 * _settings->GetRefinementThreshold() / ( _a * pow( _TRef, 4 ) ) );
    _settings->SetCoarsenThreshold( 1e-15 * _settings->GetCoarsenThreshold() / ( _a * pow( _TRef, 4 ) ) );

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ThermalRadiativeGeneral] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ThermalRadiativeGeneral::~ThermalRadiativeGeneral() {}

Vector ThermalRadiativeGeneral::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    // Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
    Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * _AbsA * ( v - u );
    g[2]     = 0.0;    // set temperature flux to zero
    return g;
}

Matrix ThermalRadiativeGeneral::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ThermalRadiativeGeneral::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = 1.0 / _c / _epsilon * u[1];
    flux( 1, 0 ) = _c / _epsilon / 3.0 * u[0];
    flux( 2, 0 ) = 0.0;    // set temperature flux to zero
    return flux;
}

Matrix ThermalRadiativeGeneral::Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const {
    unsigned nStates             = static_cast<unsigned>( uQ.rows() );
    unsigned Nq                  = _settings->GetNqPEAtRef( level );
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( level );    // get indices in quadrature array for current refinement level

    Matrix y( nStates, Nq, 0.0 );
    double S           = 0.0;    // source, needs to be defined
    double varianceVal = 0;

    for( unsigned k = 0; k < qIndex.size(); ++k ) {
        if( _suOlson && t < 10 && std::fabs( x[0] ) < 0.5 + _variances[0] * _xiQuad[qIndex[k]][0] ) {
            S           = _a;
            varianceVal = _variances[0];
        }
        else {
            S = 0.0;
        }

        double Q = S / _sigma / _a / std::pow( _TRef, 4 );

        double E      = uQ( 0, k );
        double F      = uQ( 1, k );
        double eTilde = uQ( 2, k );    // scaled internal energy
        if( eTilde < 0 ) {
            std::cout << "eTilde < 0 !!!!" << std::endl;
            std::cout << "eTilde = " << eTilde << std::endl;
            std::cout << "E = " << E << std::endl;
            std::cout << "F = " << F << std::endl;
        }
        double TTilde = ScaledTemperature( eTilde );
        // y( 0, k ) = ( -( E - U ) + ( Q + varianceVal * _xiQuad[k][0] ) ) / _epsilon;

        y( 0, k ) = ( -( E - std::pow( TTilde, 4 ) ) + Q ) / _epsilon;
        y( 1, k ) = -F / _epsilon;
        y( 2, k ) = ( E - std::pow( TTilde, 4 ) ) / _epsilon;
    }

    return y;
}

double ThermalRadiativeGeneral::ScaledInternalEnergy( double TTilde ) const {
    double T = TTilde * _TRef;
    // std::cout << "T = " << T << std::endl;
    double e;
    if( _constitutiveLaw == 1 ) {
        e = _alpha / 4.0 * pow( T, 4 );
    }
    else {
        e = _cV * T;
    }
    // std::cout << "cV * T = " << e << std::endl;
    return e / ( _a * pow( _TRef, 4 ) );
}

double ThermalRadiativeGeneral::ScaledTemperature( double eTilde ) const {
    double e = eTilde * _a * pow( _TRef, 4 );
    double T;
    if( _constitutiveLaw == 1 ) {
        T = pow( 4.0 * e / _alpha, 1.0 / 4.0 );
    }
    else {
        T = e / _cV;
    }
    return T / _TRef;
}

Matrix ThermalRadiativeGeneral::F( const Matrix& u ) {
    _log->error( "[ThermalRadiativeGeneral] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ThermalRadiativeGeneral::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    double cfl = _settings->GetCFL();

    double maxVelocity = std::sqrt( 1 / 3.0 ) / _epsilon;

    return ( cfl * dx ) / maxVelocity;
}

Vector ThermalRadiativeGeneral::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    auto sigma = _settings->GetSigma();
    double E, F, internalEnergy;
    double T = 0;

    if( _suOlson ) {
        E = 0.0;    // std::fmax( 1e-4 * _a,
                    //_a * pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) *
                    //  exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );
        F              = 0;
        internalEnergy = 1e-7 * _a * pow( _TRef, 4 );    // fix to ensure positive values of the inner energy - use 1e-3 without IPM
    }
    else {
        double a = 0.275;
        double b = 0.1;
        if( xi.size() > 1 ) {
            a = a + sigma[1] * xi[1];
        }
        if( xi.size() > 2 ) {
            b = b + sigma[2] * xi[2];
        }
        double alpha = pow( a, 1.0 / 4.0 );
        double beta  = pow( b, 1.0 / 4.0 );
        double tau0  = 0.1 + sigma[0] * xi[0];

        if( x[0] < 0.0 )
            T = alpha;
        else if( x[0] < tau0 )
            T = 1;
        else
            T = beta;
        F              = 0.0;
        internalEnergy = ScaledInternalEnergy( T / _TRef ) * ( _a * pow( _TRef, 4 ) );    // 1e-7 * _a * pow( T, 4 );
        E              = std::pow( T / _TRef, 4 );
    }

    y[0] = E / _a / pow( _TRef, 4 );
    if( !_suOlson ) y[0] = std::pow( T / _TRef, 4 );
    y[1] = F / _a / pow( _TRef, 4 );
    y[2] = internalEnergy / ( _a * pow( _TRef, 4 ) );
    // std::cout << y[2] << std::endl;
    return y;
}

Vector ThermalRadiativeGeneral::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalRadiativeGeneral: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
