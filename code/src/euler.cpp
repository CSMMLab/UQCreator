#include "euler.h"

Euler::Euler( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Euler] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Euler::~Euler() {}

void Euler::Solve() {}

Vector Euler::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    double rhoInv = 1.0 / u[0];    //
    double vU     = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[0];
    double vV     = v[1] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[2] - 0.5 * v[0] * pow( vV, 2 ) );
    double aV     = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * vU;
    double uVProjected = nUnit[0] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        return F( u ) * n;
    else if( lambdaMax <= 0 )
        return F( v ) * n;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
}

Matrix Euler::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix Euler::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v + p;
    flux( 2, 0 ) = ( u[2] + p ) * v;
    return flux;
}

Matrix Euler::F( const Matrix& u ) {
    _log->error( "[euler] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double Euler::ComputeDt( const Matrix& u, double dx ) const {
    double dtMinTotal = 1e10;
    double dtMin;
    double rhoInv, v, p, a, cfl;

    cfl = _settings->GetCFL();

    for( unsigned k = 0; k < u.columns(); ++k ) {
        rhoInv = 1.0 / u( 0, k );
        v      = u( 1, k ) * rhoInv;
        p      = ( _gamma - 1.0 ) * ( u( 2, k ) - 0.5 * u( 0, k ) * pow( v, 2 ) );
        a      = sqrt( _gamma * p * rhoInv );

        dtMin      = ( cfl / dx ) * std::min( std::fabs( 1.0 / ( v - a ) ), std::fabs( 1.0 / ( v + a ) ) );
        dtMinTotal = std::min( dtMin, dtMinTotal );
    }

    return dtMinTotal;
}

Vector Euler::IC( const Vector& x, const Vector& xi ) {
    double x0    = 0.3;
    double gamma = 1.4;

    double rhoL = 1.0;
    double rhoR = 0.125;
    double pL   = 1.0;
    double pR   = 0.1;
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
        y[0]                  = rhoR;
        y[1]                  = rhoR * uR;
        double kineticEnergyR = 0.5 * rhoR * pow( uR, 2 );
        double innerEnergyR   = ( pR / ( rhoR * ( gamma - 1 ) ) ) * rhoR;
        y[2]                  = kineticEnergyR + innerEnergyR;
    }
    return y;
}

Vector Euler::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[Euler: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
