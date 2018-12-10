#include "shallowwaterclosure.h"

ShallowWaterClosure::ShallowWaterClosure( Settings* settings ) : Closure( settings ), _g( 9.81 ) {
    _nStates = 2;
    _alpha   = 1.0;
}

ShallowWaterClosure::~ShallowWaterClosure() {}

void ShallowWaterClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = ( 0.5 * ( 2 * Lambda[0] + pow( Lambda[1], 2 ) ) ) / ( _g );
    out[1] = ( 0.5 * ( 2 * Lambda[0] * Lambda[1] + pow( Lambda[1], 3 ) ) ) / _g;
}

void ShallowWaterClosure::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        out( 0, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) + pow( Lambda( 1, k ), 2 ) ) ) / ( _g );
        out( 1, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 1, k ) + pow( Lambda( 1, k ), 3 ) ) ) / _g;
    }
}

Matrix ShallowWaterClosure::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        y( 0, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) + pow( Lambda( 1, k ), 2 ) ) ) / _g;
        y( 1, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 1, k ) + pow( Lambda( 1, k ), 3 ) ) ) / _g;
    }

    return y;
}

void ShallowWaterClosure::DU( Matrix& y, const Vector& Lambda ) {
    y( 0, 0 ) = 1.0 / _g;
    y( 0, 1 ) = Lambda[1] / _g;
    y( 1, 0 ) = Lambda[1] / _g;
    y( 1, 1 ) = ( Lambda[0] + 1.5 * pow( Lambda[1], 2 ) ) / _g;
}

void ShallowWaterClosure::DS( Vector& ds, const Vector& u ) const {
    ds[0] = ( _g * pow( u[0], 3 ) - 0.5 * pow( u[1], 2 ) ) / pow( u[0], 2 );
    ds[1] = u[1] / u[0];
}
