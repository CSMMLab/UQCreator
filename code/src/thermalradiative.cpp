#include "thermalradiative.h"

ThermalRadiative::ThermalRadiative( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );

    // physical constants
    _c     = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a     = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef  = 1.0;                    // reference temperature
    _sigma = 1.0;                    // opacity
    _alpha = 1.0;                    // heat capacity parameter c_v = alpha T^3

    _epsilon = 4 * _a / _alpha;

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
    double dtMinTotal = 1e10;
    double dtMin;
    double rhoInv, v, p, a, cfl;
    unsigned kEnd = _settings->GetNqPEAtRef( level );

    cfl = _settings->GetCFL();

    for( unsigned k = 0; k < kEnd; ++k ) {
        rhoInv = 1.0 / u( 0, k );
        v      = u( 1, k ) * rhoInv;
        p      = ( _gamma - 1.0 ) * ( u( 2, k ) - 0.5 * u( 0, k ) * pow( v, 2 ) );
        a      = sqrt( _gamma * p * rhoInv );

        dtMin      = ( cfl * dx ) * std::min( std::fabs( 1.0 / ( v - a ) ), std::fabs( 1.0 / ( v + a ) ) );
        dtMinTotal = std::min( dtMin, dtMinTotal );
    }

    return dtMinTotal;
}

Vector ThermalRadiative::IC( const Vector& x, const Vector& xi ) {
    double x0    = 0.3;
    double gamma = 1.4;

    double rhoL = 1.0;
    double rhoR = 0.1;
    double pL   = 1.0;
    double pR   = 0.125;
    double uL   = 0.0;
    double uR   = 0.0;
    Vector y( _nStates );
    _sigma = _settings->GetSigma();
    if( x[0] < x0 + _sigma[0] * xi[0] ) {
        y[0]                  = rhoL;
        y[1]                  = rhoL * uL;
        double kineticEnergyL = 0.5 * rhoL * pow( uL, 2 );
        double innerEnergyL   = ( pL / ( rhoL * ( gamma - 1 ) ) ) * rhoL;
        y[2]                  = kineticEnergyL + innerEnergyL;
    }
    else {
        y[0] = rhoR;
        if( xi.size() > 1 ) {
            y[0] += _sigma[1] * xi[1];
        }
        if( xi.size() > 2 ) {
            pR += _sigma[2] * xi[2];
        }
        y[1]                  = rhoR * uR;
        double kineticEnergyR = 0.5 * rhoR * pow( uR, 2 );
        double innerEnergyR   = ( pR / ( rhoR * ( gamma - 1 ) ) ) * rhoR;
        y[2]                  = kineticEnergyR + innerEnergyR;
    }
    return y;
}

Vector ThermalRadiative::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalRadiative: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
