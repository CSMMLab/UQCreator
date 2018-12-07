#include "shallowwaterclosure2d.h"

ShallowWaterClosure2D::ShallowWaterClosure2D( Settings* settings ) : Closure( settings ), _g( 9.81 ) {
    _nStates = 3;
    _alpha   = 0.1;
}

ShallowWaterClosure2D::~ShallowWaterClosure2D() {}

void ShallowWaterClosure2D::U( Vector& out, const Vector& Lambda ) {
    out[0] = ( 0.5 * ( 2 * Lambda[0] * Lambda[2] + pow( Lambda[1], 2 ) * Lambda[2] + pow( Lambda[2], 3 ) ) ) / ( _g * Lambda[2] );
    out[1] = ( 0.5 * Lambda[1] * ( 2 * Lambda[0] * Lambda[2] + pow( Lambda[1], 2 ) * Lambda[2] + pow( Lambda[2], 3 ) ) ) / ( _g * Lambda[2] );
    out[2] = ( ( 0.5 * ( 2 * Lambda[0] * Lambda[2] + pow( Lambda[1], 2 ) * Lambda[2] + pow( Lambda[2], 3 ) ) ) / _g );
}

void ShallowWaterClosure2D::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        out( 0, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) /
                      ( _g * Lambda( 2, k ) );
        out( 1, k ) = ( 0.5 * Lambda( 1, k ) *
                        ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) /
                      ( _g * Lambda( 2, k ) );
        out( 2, k ) =
            ( ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) / _g );
    }
}

Matrix ShallowWaterClosure2D::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        y( 0, k ) = ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) /
                    ( _g * Lambda( 2, k ) );
        y( 1, k ) = ( 0.5 * Lambda( 1, k ) *
                      ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) /
                    ( _g * Lambda( 2, k ) );
        y( 2, k ) = ( ( 0.5 * ( 2 * Lambda( 0, k ) * Lambda( 2, k ) + pow( Lambda( 1, k ), 2 ) * Lambda( 2, k ) + pow( Lambda( 2, k ), 3 ) ) ) / _g );
    }

    return y;
}

void ShallowWaterClosure2D::DU( Matrix& y, const Vector& Lambda ) {
    y( 0, 0 ) = 1.0 / _g;
    y( 0, 1 ) = Lambda[1] / _g;
    y( 0, 2 ) = Lambda[2] / _g;
    y( 1, 0 ) = Lambda[1] / _g;
    y( 1, 1 ) = ( Lambda[0] + 1.5 * pow( Lambda[1], 2 ) + 0.5 * pow( Lambda[2], 2 ) ) / _g;
    y( 1, 2 ) = ( Lambda[1] * Lambda[2] ) / _g;
    y( 2, 0 ) = Lambda[2] / _g;
    y( 2, 1 ) = ( Lambda[1] * Lambda[2] ) / _g;
    y( 2, 2 ) = ( 0.5 * ( 2.0 * Lambda[0] + pow( Lambda[1], 2 ) + 3.0 * pow( Lambda[2], 2 ) ) ) / _g;
}

void ShallowWaterClosure2D::DS( Vector& ds, const Vector& u ) const {
    ds[0] = ( _g * pow( u[0], 3 ) - 0.5 * pow( u[1], 2 ) - 0.5 * pow( u[2], 2 ) ) / pow( u[0], 2 );
    ds[1] = u[1] / u[0];
    ds[2] = u[2] / u[0];
}
