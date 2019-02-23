#include "eulerclosure2d.h"

EulerClosure2D::EulerClosure2D( Settings* settings ) : Closure( settings ), _gamma( -settings->GetGamma() ) { _alpha = 1.0; }

EulerClosure2D::~EulerClosure2D() {}

void EulerClosure2D::U( Vector& out, const Vector& Lambda ) {
    double v1      = Lambda[0];
    double v2      = Lambda[1];
    double v3      = Lambda[2];
    double v4      = Lambda[3];
    double expTerm = exp( 1.0 / ( 1.0 + _gamma ) * ( ( pow( v2, 2 ) + pow( v3, 2 ) - 2.0 * v1 * v4 - 2.0 * v4 * _gamma ) / ( 2.0 * v4 ) ) ) *
                     pow( -v4, 1.0 / ( 1.0 + _gamma ) );
    out[0] = expTerm;
    out[1] = -( ( v2 * expTerm ) / v4 );
    out[2] = -( ( v3 * expTerm ) / v4 );
    out[3] = -( ( expTerm * ( -pow( v2, 2 ) - pow( v3, 2 ) + 2.0 * v4 ) ) / ( 2.0 * pow( v4, 2 ) ) );
}

void EulerClosure2D::U( Matrix& out, const Matrix& Lambda ) {
    double expTerm, v1, v2, v3, v4;
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        v1      = Lambda( 0, k );
        v2      = Lambda( 1, k );
        v3      = Lambda( 2, k );
        v4      = Lambda( 3, k );
        expTerm = pow( -exp( ( ( pow( v2, 2 ) + pow( v3, 2 ) - 2.0 * v1 * v4 - 2.0 * v4 * _gamma ) / ( 2.0 * v4 ) ) ) * v4, 1.0 / ( 1.0 + _gamma ) );
        out( 0, k ) = expTerm;
        out( 1, k ) = -( ( v2 * expTerm ) / v4 );
        out( 2, k ) = -( ( v3 * expTerm ) / v4 );
        out( 3, k ) = -( ( expTerm * ( -pow( v2, 2 ) - pow( v3, 2 ) + 2.0 * v4 ) ) / ( 2.0 * pow( v4, 2 ) ) );
    }
}

Matrix EulerClosure2D::U( const Matrix& Lambda ) {
    double expTerm, v1, v2, v3, v4;
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        v1      = Lambda( 0, k );
        v2      = Lambda( 1, k );
        v3      = Lambda( 2, k );
        v4      = Lambda( 3, k );
        expTerm = pow( -exp( ( ( pow( v2, 2 ) + pow( v3, 2 ) - 2.0 * v1 * v4 - 2.0 * v4 * _gamma ) / ( 2.0 * v4 ) ) ) * v4, 1.0 / ( 1.0 + _gamma ) );

        y( 0, k ) = expTerm;
        y( 1, k ) = -( ( v2 * expTerm ) / v4 );
        y( 2, k ) = -( ( v3 * expTerm ) / v4 );
        y( 3, k ) = -( ( expTerm * ( -pow( v2, 2 ) - pow( v3, 2 ) + 2.0 * v4 ) ) / ( 2.0 * pow( v4, 2 ) ) );
    }

    return y;
}

void EulerClosure2D::DU( Matrix& y, const Vector& Lambda ) {
    double v1        = Lambda[0];
    double v2        = Lambda[1];
    double v3        = Lambda[2];
    double v4        = Lambda[3];
    double v4pow3    = pow( v4, 3 );
    double v4pow2    = pow( v4, 2 );
    double v4pow2Inv = 1.0 / v4pow2;
    double v3pow2    = pow( v3, 2 );
    double v2pow2    = pow( v2, 2 );
    // double expTerm =
    //    pow( -exp( ( ( v2pow2 + v3pow2 - 2.0 * v4 * ( v1 + _gamma ) ) / ( 2.0 * v4 ) ) ) * v4, 1.0 / ( 1.0 + _gamma ) ) * 1.0 / ( 1.0 + _gamma );

    double expTerm = exp( 1.0 / ( 1.0 + _gamma ) * ( ( pow( v2, 2 ) + pow( v3, 2 ) - 2.0 * v1 * v4 - 2.0 * v4 * _gamma ) / ( 2.0 * v4 ) ) ) *
                     pow( -v4, 1.0 / ( 1.0 + _gamma ) ) / ( 1.0 + _gamma );
    double vTerm = expTerm * ( v2pow2 + v3pow2 + 2.0 * v4 * _gamma ) / ( 2 * v4pow3 );
    y( 0, 0 )    = -( expTerm );
    y( 0, 1 )    = ( v2 * expTerm ) / ( v4 );
    y( 0, 2 )    = ( v3 * expTerm ) / ( v4 );
    y( 0, 3 )    = -( 0.5 * ( ( v2pow2 + v3pow2 - 2.0 * v4 ) * expTerm ) * v4pow2Inv );
    y( 1, 0 )    = y( 0, 1 );
    y( 1, 1 )    = -( ( expTerm * ( v2pow2 + ( _gamma + 1.0 ) * v4 ) ) * v4pow2Inv );
    y( 1, 2 )    = -( ( v2 * v3 * expTerm ) * v4pow2Inv );
    y( 1, 3 )    = v2 * vTerm;
    y( 2, 0 )    = y( 0, 2 );
    y( 2, 1 )    = y( 1, 2 );
    y( 2, 2 )    = -( ( expTerm * ( v3pow2 + ( _gamma + 1.0 ) * v4 ) ) * v4pow2Inv );
    y( 2, 3 )    = v3 * vTerm;
    y( 3, 0 )    = y( 0, 3 );
    y( 3, 1 )    = y( 1, 3 );
    y( 3, 2 )    = y( 2, 3 );
    y( 3, 3 )    = -( expTerm * ( pow( v2, 4 ) + pow( v3, 4 ) + 4.0 * v3pow2 * v4 * _gamma - 4.0 * v4pow2 * _gamma +
                               2.0 * v2pow2 * ( v3pow2 + 2.0 * v4 * _gamma ) ) ) /
                ( 4.0 * pow( v4, 4 ) );
}

void EulerClosure2D::DS( Vector& ds, const Vector& u ) const {
    double gamma = _gamma;
    double rho   = u[0];
    double rhoU  = u[1];
    double rhoV  = u[2];
    double rhoV2 = pow( rhoV, 2 );
    double rhoU2 = pow( rhoU, 2 );
    double rhoE  = u[3];
    ds[0]        = ( rhoU2 + rhoV2 + gamma * ( 2 * rho * rhoE - rhoU2 - rhoV2 ) ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) -
            std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );
    ds[1] = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[2] = -( ( 2 * rho * rhoV ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[3] = -( rho / ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );

    // safety factor on ds3
    // ds[3] -= 1e-7;
}
