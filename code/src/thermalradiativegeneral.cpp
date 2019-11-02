#include "thermalradiativegeneral.h"
#include "quadraturegrid.h"

ThermalRadiativeGeneral::ThermalRadiativeGeneral( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    _settings->SetSource( true );

    // physical constants
    _c             = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a             = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef          = 1.0;                    // reference temperature
    _sigma         = 1.0;                    // opacity
    _alpha         = 4.0 * _a;               // heat capacity parameter c_v = alpha T^3
    double sigmaSB = 5.6704 * 1e-5;          // Stefan Boltzmann constant in [erg/cm^2/s/K^4]
    _a             = 4.0 * sigmaSB / _c;
    //_c             = 1.0;
    //_a             = 1.0;

    _epsilon = 4.0 * _a / _alpha;

    // compute xi Quadrature points
    Vector xiEta( _settings->GetNDimXi() );
    _variances = _settings->GetSigma();

    auto grid = QuadratureGrid::Create( _settings, _settings->GetNQTotal() );
    _xiQuad   = grid->GetNodes();

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
    Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
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
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq, 0.0 );
    double S           = 0.0;    // source, needs to be defined
    double varianceVal = 0;

    // std::cout << "level " << level << ", Nq = " << Nq << std::endl;

    for( unsigned k = 0; k < Nq; ++k ) {
        if( t < 10 && std::fabs( x[0] ) < 0.5 + _variances[0] * _xiQuad[k][0] ) {
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
        y( 2, k ) = E - std::pow( TTilde, 4 );
    }

    return y;
}

double ThermalRadiativeGeneral::ScaledInternalEnergy( double TTilde ) const {
    double T = TTilde * _TRef;
    double e = _alpha / 4.0 * pow( T, 4 );
    return e / ( _a * pow( _TRef, 4 ) );
}

double ThermalRadiativeGeneral::ScaledTemperature( double eTilde ) const {
    double e = eTilde * _a * pow( _TRef, 4 );
    double T = pow( 4.0 * e / _alpha, 1.0 / 4.0 );
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
    auto sigma     = _settings->GetSigma();
    double sigmaXi = sigma[0] * xi[0];

    double E = 0.0;    // std::fmax( 1e-4 * _a,
                       //_a * pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) *
                       //  exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );
    double F              = 0;
    double internalEnergy = 1e-7 * _a * pow( _TRef, 4 );    // fix to ensure positive values of the inner energy - use 1e-3 without IPM

    y[0] = E / _a / pow( _TRef, 4 );
    y[1] = F / _a / pow( _TRef, 4 );
    y[2] = internalEnergy / ( _a * pow( _TRef, 4 ) );
    return y;
}

Vector ThermalRadiativeGeneral::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalRadiativeGeneral: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
