#include "eulerclosure.h"

EulerClosure::EulerClosure( Problem* problem ) : Closure( problem ), _gamma( problem->GetGamma() ) {}

EulerClosure::~EulerClosure() {}

void EulerClosure::U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda ) {
    out[1] = exp( ( 2.0 * Lambda[1] * Lambda[3] - 2.0 * Lambda[3] * log( -Lambda[3] ) - 2.0 * Lambda[3] * _gamma - pow( Lambda[2], 2.0 ) ) /
                  ( 2.0 * Lambda[3] * ( _gamma - 1.0 ) ) );
    out[2] = -( Lambda[2] / Lambda[3] ) * out[1];
    out[3] = ( ( pow( Lambda[2], 2.0 ) - 2.0 * Lambda[3] ) / ( 2.0 * pow( Lambda[3], 2.0 ) ) ) * out[1];
}

blaze::DynamicMatrix<double> EulerClosure::U( const blaze::DynamicMatrix<double>& Lambda ) {
    blaze::DynamicMatrix<double> y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
        y( 1, k ) = exp( ( 2.0 * Lambda( 1, k ) * Lambda( 2, k ) - 2.0 * Lambda( 3, k ) * log( -Lambda( 3, k ) ) - 2.0 * Lambda( 3, k ) * _gamma -
                           pow( Lambda( 2, k ), 2.0 ) ) /
                         ( 2.0 * Lambda( 3, k ) * ( _gamma - 1.0 ) ) );
        y( 2, k ) = -( Lambda( 2, k ) / Lambda( 3, k ) ) * y( 1, k );
        y( 3, k ) = ( ( pow( Lambda( 2, k ), 2.0 ) - 2.0 * Lambda( 3, k ) ) / ( 2.0 * pow( Lambda( 3, k ), 2.0 ) ) ) * y( 1, k );
    }

    return y;
}

void EulerClosure::DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda ) {
    double E     = exp( ( 2.0 * Lambda[1] * Lambda[3] - 2.0 * Lambda[3] * log( -Lambda[3] ) - 2.0 * Lambda[3] * _gamma - pow( Lambda[2], 2 ) ) /
                    ( 2.0 * Lambda[3] * ( _gamma - 1.0 ) ) );
    double dEdv1 = ( 1.0 / ( _gamma - 1.0 ) ) * E;
    double dEdv2 = -( Lambda[2] / ( Lambda[3] * ( _gamma - 1.0 ) ) ) * E;
    double dEdv3 = ( ( -2.0 * _gamma + 2.0 * Lambda[1] - 2.0 * log( -Lambda[3] ) - 2.0 ) / ( 2.0 * ( _gamma - 1.0 ) * Lambda[3] ) -
                     ( -2.0 * _gamma * Lambda[3] + 2.0 * Lambda[1] * Lambda[3] - pow( Lambda[2], 2 ) - 2.0 * Lambda[3] * log( -Lambda[3] ) ) /
                         ( 2.0 * ( _gamma - 1.0 ) * pow( Lambda[3], 2 ) ) ) *
                   E;
    y( 1, 1 ) = dEdv1;
    y( 1, 2 ) = dEdv2;
    y( 1, 3 ) = dEdv3;
    y( 2, 1 ) = -( Lambda[2] / Lambda[3] ) * dEdv1;
    y( 2, 2 ) = -( 1.0 / Lambda[3] ) * E - ( Lambda[2] / Lambda[3] ) * dEdv2;
    y( 2, 3 ) = ( Lambda[2] / ( pow( Lambda[3], 2 ) ) ) * E - ( Lambda[2] / Lambda[3] ) * dEdv3;
    y( 3, 1 ) = ( ( pow( Lambda[2], 2 ) - 2.0 * Lambda[3] ) / ( 2.0 * pow( Lambda[3], 2 ) ) ) * dEdv1;
    y( 3, 2 ) = ( ( pow( Lambda[2], 2 ) - 2.0 * Lambda[3] ) / ( 2.0 * pow( Lambda[3], 2 ) ) ) * dEdv2 +
                ( 2.0 * Lambda[2] / ( 2.0 * pow( Lambda[3], 2 ) ) ) * E;
    y( 3, 3 ) = ( pow( Lambda[2], 2 ) / ( pow( Lambda[3], -3 ) ) + pow( Lambda[3], -2 ) ) * E +
                ( ( pow( Lambda[2], 2.0 ) - 2.0 * Lambda[3] ) / ( 2.0 * pow( Lambda[3], 2 ) ) ) * dEdv3;
}
