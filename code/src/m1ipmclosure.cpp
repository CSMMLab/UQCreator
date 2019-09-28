#include "m1ipmclosure.h"

M1IPMClosure::M1IPMClosure( Settings* settings ) : Closure( settings ), _gamma( _settings->GetGamma() ) { _alpha = 1.0; }

M1IPMClosure::~M1IPMClosure() {}

void M1IPMClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = ( -exp( Lambda[0] - Lambda[1] ) + exp( Lambda[0] + Lambda[1] ) ) / Lambda[1];
    out[1] = ( exp( Lambda[0] + Lambda[1] ) * ( -1.0 + Lambda[1] ) + exp( Lambda[0] - Lambda[1] ) * ( 1.0 + Lambda[1] ) ) / pow( Lambda[1], 2 );
}

void M1IPMClosure::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        out( 0, k ) = ( -exp( Lambda( 0, k ) - Lambda( 1, k ) ) + exp( Lambda( 0, k ) + Lambda( 1, k ) ) ) / Lambda( 1, k );
        out( 1, k ) = ( exp( Lambda( 0, k ) + Lambda( 1, k ) ) * ( -1.0 + Lambda( 1, k ) ) +
                        exp( Lambda( 0, k ) - Lambda( 1, k ) ) * ( 1.0 + Lambda( 1, k ) ) ) /
                      pow( Lambda( 1, k ), 2 );
    }
}

Matrix M1IPMClosure::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        y( 0, k ) = ( -exp( Lambda( 0, k ) - Lambda( 1, k ) ) + exp( Lambda( 0, k ) + Lambda( 1, k ) ) ) / Lambda( 1, k );
        y( 1, k ) = ( exp( Lambda( 0, k ) + Lambda( 1, k ) ) * ( -1.0 + Lambda( 1, k ) ) +
                      exp( Lambda( 0, k ) - Lambda( 1, k ) ) * ( 1.0 + Lambda( 1, k ) ) ) /
                    pow( Lambda( 1, k ), 2 );
    }

    return y;
}

void M1IPMClosure::DU( Matrix& y, const Vector& Lambda ) {
    y( 0, 0 ) = ( exp( Lambda[0] - Lambda[1] ) * ( -1.0 + exp( 2.0 * Lambda[1] ) ) ) / Lambda[1];
    y( 0, 1 ) = ( exp( Lambda[0] - Lambda[1] ) * ( 1.0 + exp( 2.0 * Lambda[1] ) * ( -1.0 + Lambda[1] ) + Lambda[1] ) ) / pow( Lambda[1], 2 );
    y( 1, 0 ) = y( 0, 1 );
    y( 1, 1 ) = ( exp( Lambda[0] - Lambda[1] ) *
                  ( -2.0 - 2.0 * Lambda[1] - pow( Lambda[1], 2 ) + exp( 2.0 * Lambda[1] ) * ( 2.0 - 2.0 * Lambda[1] + pow( Lambda[1], 2 ) ) ) ) /
                pow( Lambda[1], 3 );
    std::cout << y << std::endl;
}

void M1IPMClosure::DS( Vector& ds, const Vector& u ) const {
    ds[0] = 5.0;
    ds[1] = 0.01;
}
