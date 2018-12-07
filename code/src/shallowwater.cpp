#include "shallowwater.h"

ShallowWater::ShallowWater( Settings* settings ) : Problem( settings ) {
    _nStates = 2;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _g           = 9.81;
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ShallowWater] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ShallowWater::~ShallowWater() {}

void ShallowWater::Solve() {}

Vector ShallowWater::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    /*double hAvg  = 0.5 * ( u[0] + v[0] );
    double h2Avg = 0.5 * ( pow( u[0], 2 ) + pow( v[0], 2 ) );
    double huAvg = 0.5 * ( u[1] + v[1] );
    double uU    = u[1] / u[0];
    double vU    = v[1] / v[0];
    double uAvg  = 0.5 * ( uU + vU );

    // Vector y( 2 );

    // y[0] = huAvg;
    // y[1] = huAvg*uAvg + _g * hAvg - 0.5*_g*h2Avg;

    Matrix y = 0.5 * ( F( u ) + F( v ) );

    // return y * n[0];
    return y * n;*/
    double uU = u[1] / u[0];
    double uV = v[1] / v[0];
    double cU = sqrt( _g * u[0] );
    double cV = sqrt( _g * v[0] );

    double lambdaMin = uU - cU;
    double lambdaMax = uV + cV;

    if( lambdaMin >= 0 )
        return F( u ) * n;
    else if( lambdaMax <= 0 )
        return F( v ) * n;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
}

Matrix ShallowWater::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ShallowWater::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = pow( u[1], 2 ) / u[0] + 0.5 * _g * pow( u[0], 2 );
    return flux;
}

Matrix ShallowWater::F( const Matrix& u ) {
    _log->error( "[ShallowWater] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ShallowWater::ComputeDt( Vector& u, double dx ) const {
    _log->error( "[ShallowWater] ComputeDt not implemented" );
    return 0.0;
}
