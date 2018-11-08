#include "burgers.h"

Burgers::Burgers( Settings* settings ) : Problem( settings ) {
    _nStates = 1;
    _settings->SetNStates( _nStates );
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[burgers] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Vector Burgers::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    if( u[0] * nUnit[0] > 0 ) {
        return F( u[0] ) * n[0];
    }
    else {
        return F( v[0] ) * n[0];
    }
}

Matrix Burgers::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = u.rows();
    unsigned Nq      = u.columns();
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Vector Burgers::F( double u ) {
    Vector y( 1 );
    y[0] = 0.5 * u * u;
    return y;
}

Matrix Burgers::F( const Matrix& u ) { return 0.5 * pow( u, 2 ); }

double Burgers::IC( double x, double uL, double uR ) {
    double a = 0.5;
    double b = 1.5;
    if( x < a ) {
        return uL;
    }
    else if( x > a && x < b ) {
        return uR + ( uL - uR ) * ( b - x ) / ( b - a );
    }
    else {
        return uR;
    }
}

double Burgers::ComputeDt( Vector& u, double dx ) const { return dx * _settings->GetCFL() / u[0]; }
