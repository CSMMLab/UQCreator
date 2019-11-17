#include "thermalradiative.h"
#include "quadraturegrid.h"

ThermalRadiative::ThermalRadiative( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    _settings->SetSource( true );
    _suOlson = false;

    // physical constants
    _c             = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a             = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef          = 1.0;                    // reference temperature
    _sigma         = 1.0;                    // opacity (controls interaction strength between material and particles) in the general case
    _alpha         = 4.0 * _a;               // heat capacity parameter c_v = alpha T^3
    double sigmaSB = 5.6704 * 1e-5;          // Stefan Boltzmann constant in [erg/cm^2/s/K^4]
    _a             = 4.0 * sigmaSB / _c;
    //_c             = 1.0;
    //_a             = 1.0;

    _epsilon = 4.0 * _a / _alpha;    // closure relation c_v(T) = alpha*T^3

    // compute xi Quadrature points
    Vector xiEta( _settings->GetNDimXi() );
    _variances = _settings->GetSigma();

    // get quadrature grid
    auto grid = QuadratureGrid::Create( _settings, _settings->GetNQTotal() );
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

    // std::cout << PInv * LambdaAbs * P << std::endl;
    // std::cout << 1.0 / _c / _epsilon << std::endl;
    // std::cout << _c / _epsilon / 3.0 << std::endl;
    // exit( EXIT_FAILURE );

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ThermalRadiative] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ThermalRadiative::~ThermalRadiative() {}

Vector ThermalRadiative::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    // Lax-Friedrichs
    // Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
    // upwinding
    Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * _AbsA * ( v - u );
    g[2]     = 0.0;    // set temperature flux to zero
    return g;
}

Matrix ThermalRadiative::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ThermalRadiative::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = 1.0 / _c / _epsilon * u[1];
    flux( 1, 0 ) = _c / _epsilon / 3.0 * u[0];
    flux( 2, 0 ) = 0.0;    // set temperature flux to zero
    return flux;
}

Matrix ThermalRadiative::Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const {
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq, 0.0 );
    double S           = 0.0;    // source, needs to be defined
    double varianceVal = 0;

    // std::cout << "level " << level << ", Nq = " << Nq << std::endl;

    for( unsigned k = 0; k < Nq; ++k ) {
        if( _suOlson && t < 10 && std::fabs( x[0] ) < 0.5 + _variances[0] * _xiQuad[k][0] ) {
            S           = _a;
            varianceVal = _variances[0];
        }
        else {
            S = 0.0;
        }

        double Q = S / _sigma / _a / std::pow( _TRef, 4 );

        double E = uQ( 0, k );
        double F = uQ( 1, k );
        double U = uQ( 2, k );
        // y( 0, k ) = ( -( E - U ) + ( Q + varianceVal * _xiQuad[k][0] ) ) / _epsilon;
        y( 0, k ) = ( -( E - U ) + Q ) / _epsilon;
        y( 1, k ) = -F / _epsilon;
        y( 2, k ) = E - U;
    }

    return y;
}

Matrix ThermalRadiative::F( const Matrix& u ) {
    _log->error( "[ThermalRadiative] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ThermalRadiative::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    double cfl = _settings->GetCFL();

    double maxVelocity = std::sqrt( 1.0 / 3.0 ) / _epsilon;

    return ( cfl * dx ) / maxVelocity;
}

Vector ThermalRadiative::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    auto sigma     = _settings->GetSigma();
    double sigmaXi = sigma[0] * xi[0];
    double E, F, T;
    F = 0;

    if( _suOlson ) {
        E = 0.0;
        T = 0.0;
    }
    else {
        double a     = 0.271;
        double b     = 0.1;
        double alpha = pow( a, 1.0 / 4.0 );
        double beta  = pow( b, 1.0 / 4.0 );
        double tau0  = 0.05;
        if( x[0] < 0.0 )
            T = alpha * _TRef;
        else if( x[0] < tau0 )
            T = _TRef;
        else
            T = beta * _TRef;
        E = _a * std::pow( T, 4 );
    }

    y[0] = E / _a / pow( _TRef, 4 );
    y[1] = F / _a / pow( _TRef, 4 );
    y[2] = pow( T, 4 ) / pow( _TRef, 4 );
    return y;
}

Vector ThermalRadiative::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalRadiative: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
