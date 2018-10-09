#include "eulerclosure.h"

EulerClosure::EulerClosure( Problem* problem ) : Closure( problem ), _gamma( problem->GetGamma() ) {}

EulerClosure::~EulerClosure() {}

void EulerClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = exp( ( 2.0 * Lambda[0] * Lambda[2] - 2.0 * Lambda[2] * log( -Lambda[2] ) - 2.0 * Lambda[2] * _gamma - pow( Lambda[1], 2.0 ) ) /
                  ( 2.0 * Lambda[2] * ( _gamma - 1.0 ) ) );
    out[1] = -( Lambda[1] / Lambda[2] ) * out[0];
    out[2] = ( ( pow( Lambda[1], 2.0 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2.0 ) ) ) * out[0];
}

Matrix EulerClosure::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        y( 0, k ) = exp( ( 2.0 * Lambda( 0, k ) * Lambda( 2, k ) - 2.0 * Lambda( 2, k ) * log( -Lambda( 2, k ) ) - 2.0 * Lambda( 2, k ) * _gamma -
                           pow( Lambda( 1, k ), 2.0 ) ) /
                         ( 2.0 * Lambda( 2, k ) * ( _gamma - 1.0 ) ) );
        y( 1, k ) = -( Lambda( 1, k ) / Lambda( 2, k ) ) * y( 0, k );
        y( 2, k ) = ( ( pow( Lambda( 1, k ), 2.0 ) - 2.0 * Lambda( 2, k ) ) / ( 2.0 * pow( Lambda( 2, k ), 2.0 ) ) ) * y( 0, k );
    }

    return y;
}

void EulerClosure::DU( Matrix& y, const Vector& Lambda ) {
    double E     = exp( ( 2.0 * Lambda[0] * Lambda[2] - 2.0 * Lambda[2] * log( -Lambda[2] ) - 2.0 * Lambda[2] * _gamma - pow( Lambda[1], 2 ) ) /
                    ( 2.0 * Lambda[2] * ( _gamma - 1.0 ) ) );
    double dEdv1 = ( 1.0 / ( _gamma - 1.0 ) ) * E;
    double dEdv2 = -( Lambda[1] / ( Lambda[2] * ( _gamma - 1.0 ) ) ) * E;
    double dEdv3 = ( ( -2.0 * _gamma + 2.0 * Lambda[0] - 2.0 * log( -Lambda[2] ) - 2.0 ) / ( 2.0 * ( _gamma - 1.0 ) * Lambda[2] ) -
                     ( -2.0 * _gamma * Lambda[2] + 2.0 * Lambda[0] * Lambda[2] - pow( Lambda[1], 2 ) - 2.0 * Lambda[2] * log( -Lambda[2] ) ) /
                         ( 2.0 * ( _gamma - 1.0 ) * pow( Lambda[2], 2 ) ) ) *
                   E;
    y( 0, 0 ) = dEdv1;
    y( 0, 1 ) = dEdv2;
    y( 0, 2 ) = dEdv3;
    y( 1, 0 ) = -( Lambda[1] / Lambda[2] ) * dEdv1;
    y( 1, 1 ) = -( 1.0 / Lambda[2] ) * E - ( Lambda[1] / Lambda[2] ) * dEdv2;
    y( 1, 2 ) = ( Lambda[1] / ( pow( Lambda[2], 2 ) ) ) * E - ( Lambda[1] / Lambda[2] ) * dEdv3;
    y( 2, 0 ) = ( ( pow( Lambda[1], 2 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2 ) ) ) * dEdv1;
    y( 2, 1 ) = ( ( pow( Lambda[1], 2 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2 ) ) ) * dEdv2 +
                ( 2.0 * Lambda[1] / ( 2.0 * pow( Lambda[2], 2 ) ) ) * E;
    y( 2, 2 ) = ( pow( Lambda[1], 2 ) / ( pow( Lambda[2], -3 ) ) + pow( Lambda[2], -2 ) ) * E +
                ( ( pow( Lambda[1], 2 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2 ) ) ) * dEdv3;
}
