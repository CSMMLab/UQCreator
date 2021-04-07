#include "closures/shallowwaterclosure2d.h"

ShallowWaterClosure2D::ShallowWaterClosure2D( Settings* settings ) : Closure( settings ), _g( 9.81 ) {
    _nStates = 3;
    _alpha   = 1.0;
}

ShallowWaterClosure2D::~ShallowWaterClosure2D() {}

void ShallowWaterClosure2D::U( Vector& out, const Vector& Lambda ) {
    out[0] = ( 0.5 * ( 2 * Lambda[0] + pow( Lambda[1], 2 ) + pow( Lambda[2], 2 ) ) ) / ( _g );
    out[1] = ( 0.5 * Lambda[1] * ( 2 * Lambda[0] + pow( Lambda[1], 2 ) + pow( Lambda[2], 2 ) ) ) / ( _g );
    out[2] = ( ( 0.5 * ( 2 * Lambda[0] * Lambda[2] + pow( Lambda[1], 2 ) * Lambda[2] + pow( Lambda[2], 3 ) ) ) / _g );
}

void ShallowWaterClosure2D::U( Tensor& out, const Tensor& Lambda ) {
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            out( 0, l, k ) = ( 0.5 * ( 2 * Lambda( 0, l, k ) + pow( Lambda( 1, l, k ), 2 ) + pow( Lambda( 2, l, k ), 2 ) ) ) / ( _g );
            out( 1, l, k ) =
                ( 0.5 * Lambda( 1, l, k ) * ( 2 * Lambda( 0, l, k ) + pow( Lambda( 1, l, k ), 2 ) + pow( Lambda( 2, l, k ), 2 ) ) ) / ( _g );
            out( 2, l, k ) = ( ( 0.5 * ( 2 * Lambda( 0, l, k ) * Lambda( 2, l, k ) + pow( Lambda( 1, l, k ), 2 ) * Lambda( 2, l, k ) +
                                         pow( Lambda( 2, l, k ), 3 ) ) ) /
                               _g );
        }
    }
}

Tensor ShallowWaterClosure2D::U( const Tensor& Lambda ) {
    Tensor y( _nStates, _nMultiElements, Lambda.columns(), 0.0 );
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            y( 0, l, k ) = ( 0.5 * ( 2 * Lambda( 0, l, k ) + pow( Lambda( 1, l, k ), 2 ) + pow( Lambda( 2, l, k ), 2 ) ) ) / ( _g );
            y( 1, l, k ) =
                ( 0.5 * Lambda( 1, l, k ) * ( 2 * Lambda( 0, l, k ) + pow( Lambda( 1, l, k ), 2 ) + pow( Lambda( 2, l, k ), 2 ) ) ) / ( _g );
            y( 2, l, k ) = ( ( 0.5 * ( 2 * Lambda( 0, l, k ) * Lambda( 2, l, k ) + pow( Lambda( 1, l, k ), 2 ) * Lambda( 2, l, k ) +
                                       pow( Lambda( 2, l, k ), 3 ) ) ) /
                             _g );
        }
    }

    return y;
}

void ShallowWaterClosure2D::DU( Matrix& y, const Vector& Lambda ) {
    y( 0, 0 ) = 1.0 / _g;
    y( 0, 1 ) = Lambda[1] / _g;
    y( 0, 2 ) = Lambda[2] / _g;
    y( 1, 0 ) = y( 0, 1 );
    y( 1, 1 ) = ( Lambda[0] + 1.5 * pow( Lambda[1], 2 ) + 0.5 * pow( Lambda[2], 2 ) ) / _g;
    y( 1, 2 ) = ( Lambda[1] * Lambda[2] ) / _g;
    y( 2, 0 ) = y( 0, 2 );
    y( 2, 1 ) = y( 1, 2 );
    y( 2, 2 ) = ( 0.5 * ( 2.0 * Lambda[0] + pow( Lambda[1], 2 ) + 3.0 * pow( Lambda[2], 2 ) ) ) / _g;
}

void ShallowWaterClosure2D::DS( Vector& ds, const Vector& u ) const {
    ds[0] = ( _g * pow( u[0], 3 ) - 0.5 * pow( u[1], 2 ) - 0.5 * pow( u[2], 2 ) ) / pow( u[0], 2 );
    ds[1] = u[1] / u[0];
    ds[2] = u[2] / u[0];
}
