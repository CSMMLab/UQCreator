#include "euler2d.h"

Euler2D::Euler2D( Settings* settings ) : Problem( settings ) {
    _nStates = 4;
    _settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[euler2d] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Euler2D::~Euler2D() {}

void Euler2D::Solve() {}

Vector Euler2D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    double rhoInv = 1.0 / u[0];
    double uU     = u[1] * rhoInv;
    double vU     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
    double aU     = sqrt( _gamma * p * rhoInv );

    rhoInv    = 1.0 / v[0];
    double uV = v[1] * rhoInv;
    double vV = v[2] * rhoInv;
    p         = ( _gamma - 1.0 ) * ( v[3] - 0.5 * v[0] * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
    double aV = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
    double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

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

Matrix Euler2D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        // column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );

        double rhoInv = 1.0 / u( 0, k );
        double uU     = u( 1, k ) * rhoInv;
        double vU     = u( 2, k ) * rhoInv;
        double p      = ( _gamma - 1.0 ) * ( u( 3, k ) - 0.5 * u( 0, k ) * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
        double aU     = sqrt( _gamma * p * rhoInv );

        rhoInv    = 1.0 / v( 0, k );
        double uV = v( 1, k ) * rhoInv;
        double vV = v( 2, k ) * rhoInv;
        p         = ( _gamma - 1.0 ) * ( v( 3, k ) - 0.5 * v( 0, k ) * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
        double aV = sqrt( _gamma * p * rhoInv );

        double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
        double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

        double lambdaMin = uUProjected - aU;
        double lambdaMax = uVProjected + aV;

        if( lambdaMin >= 0 )
            column( y, k ) = F( column( u, k ) ) * n;
        else if( lambdaMax <= 0 )
            column( y, k ) = F( column( v, k ) ) * n;
        else {
            column( y, k ) = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * F( column( u, k ) ) - lambdaMin * F( column( v, k ) ) ) * n +
                                                                     lambdaMax * lambdaMin * ( column( v, k ) - column( u, k ) ) * norm( n ) );
        }
    }
    return y;
}

Matrix Euler2D::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v1     = u[1] * rhoInv;
    double v2     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( v1, 2 ) + pow( v2, 2 ) ) );
    Matrix flux( u.size(), 2 );

    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v1 + p;
    flux( 2, 0 ) = u[1] * v2;
    flux( 3, 0 ) = ( u[3] + p ) * v1;
    flux( 0, 1 ) = u[2];
    flux( 1, 1 ) = u[2] * v1;
    flux( 2, 1 ) = u[2] * v2 + p;
    flux( 3, 1 ) = ( u[3] + p ) * v2;

    return flux;
}

Matrix Euler2D::F( const Matrix& u ) {
    _log->error( "[euler2d] Flux not implemented" );
    exit( EXIT_FAILURE );
    return 0.5 * pow( u, 2 );
}

double Euler2D::ComputeDt( Vector& u, double dx ) const {
    double rhoInv = 1.0 / u[0];
    double uU     = u[1] * rhoInv;
    double vU     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
    double a      = sqrt( _gamma * p * rhoInv );

    double dt1 = dx * _settings->GetCFL() * ( uU - a );
    double dt2 = dx * _settings->GetCFL() * ( uU + a );
    double dt3 = dx * _settings->GetCFL() * ( vU - a );
    double dt4 = dx * _settings->GetCFL() * ( vU + a );
    if( dt1 < dt2 && dt1 < dt3 && dt1 < dt4 )
        return dt1;
    else if( dt2 < dt1 && dt2 < dt3 && dt2 < dt4 )
        return dt2;
    else if( dt3 < dt1 && dt3 < dt2 && dt3 < dt4 )
        return dt3;
    else
        return dt4;
}
