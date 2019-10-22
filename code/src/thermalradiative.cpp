#include "thermalradiative.h"

ThermalRadiative::ThermalRadiative( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    _settings->SetSource( false );

    // physical constants
    _c     = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a     = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef  = 1.0;                    // reference temperature
    _sigma = 1.0;                    // opacity
    _alpha = 1.0;                    // heat capacity parameter c_v = alpha T^3

    _epsilon = 4.0 * _a / _alpha;

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
    Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
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

Matrix ThermalRadiative::Source( const Matrix& uQ ) const {
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = static_cast<unsigned>( uQ.columns() );
    Matrix y( nStates, Nq, 0.0 );
    double Q = 0.0;    // source, needs to be defined

    for( unsigned k = 0; k < Nq; ++k ) {
        double E  = uQ( 0, k );
        double F  = uQ( 1, k );
        double U  = uQ( 2, k );
        y( 0, k ) = ( -( E - U ) + Q ) / _epsilon;
        y( 1, k ) = -F / _epsilon;
        y( 0, k ) = E - U;
    }

    return y;
}

Matrix ThermalRadiative::F( const Matrix& u ) {
    _log->error( "[ThermalRadiative] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ThermalRadiative::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    double cfl = _settings->GetCFL();

    double maxVelocity = 1.0 / _c / _epsilon;

    return ( cfl * dx ) / maxVelocity;
}

Vector ThermalRadiative::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double x0      = 0.0;
    double floor   = 1e-4;
    auto sigma     = _settings->GetSigma();
    double sigmaXi = sigma[0] * xi[0];

    double E = std::fmax(
        floor, pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) * exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );
    double F = 0;
    double T = 1.0;

    y[0] = E / _a / pow( _TRef, 4 );
    y[1] = F / _a / pow( _TRef, 4 );
    y[2] = pow( T, 4 ) / pow( _TRef, 4 );
    return y;
}

Vector ThermalRadiative::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalRadiative: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
