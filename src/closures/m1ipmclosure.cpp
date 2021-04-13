#include "closures/m1ipmclosure.h"

M1IPMClosure::M1IPMClosure( Settings* settings ) : Closure( settings ), _gamma( _settings->GetGamma() ) { _alpha = 1.0; }

M1IPMClosure::~M1IPMClosure() {}

void M1IPMClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = ( -exp( Lambda[0] - Lambda[1] ) + exp( Lambda[0] + Lambda[1] ) ) / Lambda[1];
    out[1] = ( exp( Lambda[0] + Lambda[1] ) * ( -1.0 + Lambda[1] ) + exp( Lambda[0] - Lambda[1] ) * ( 1.0 + Lambda[1] ) ) / pow( Lambda[1], 2 );
}

void M1IPMClosure::U( Tensor& out, const Tensor& Lambda ) {
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            out( 0, l, k ) = ( -exp( Lambda( 0, l, k ) - Lambda( 1, l, k ) ) + exp( Lambda( 0, l, k ) + Lambda( 1, l, k ) ) ) / Lambda( 1, l, k );
            out( 1, l, k ) = ( exp( Lambda( 0, l, k ) + Lambda( 1, l, k ) ) * ( -1.0 + Lambda( 1, l, k ) ) +
                               exp( Lambda( 0, l, k ) - Lambda( 1, l, k ) ) * ( 1.0 + Lambda( 1, l, k ) ) ) /
                             pow( Lambda( 1, l, k ), 2 );
        }
    }
}

Tensor M1IPMClosure::U( const Tensor& Lambda ) {
    Tensor y( _nStates, _nMultiElements, Lambda.columns(), 0.0 );
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            y( 0, l, k ) = ( -exp( Lambda( 0, l, k ) - Lambda( 1, l, k ) ) + exp( Lambda( 0, l, k ) + Lambda( 1, l, k ) ) ) / Lambda( 1, l, k );
            y( 1, l, k ) = ( exp( Lambda( 0, l, k ) + Lambda( 1, l, k ) ) * ( -1.0 + Lambda( 1, l, k ) ) +
                             exp( Lambda( 0, l, k ) - Lambda( 1, l, k ) ) * ( 1.0 + Lambda( 1, l, k ) ) ) /
                           pow( Lambda( 1, l, k ), 2 );
        }
    }

    return y;
}

void M1IPMClosure::DU( Matrix& y, const Vector& Lambda ) {
    double v0 = Lambda[0];
    double v1 = Lambda[1];

    y( 0, 0 ) = ( exp( Lambda[0] - Lambda[1] ) * ( -1.0 + exp( 2.0 * Lambda[1] ) ) ) / Lambda[1];

    y( 0, 1 ) = ( exp( v0 - v1 ) * ( 1.0 + exp( 2.0 * v1 ) * ( -1.0 + v1 ) + v1 ) ) / pow( v1, 2 );

    y( 1, 0 ) = y( 0, 1 );

    y( 1, 1 ) = ( exp( v0 - v1 ) * ( -2.0 - 2.0 * v1 - pow( v1, 2 ) + exp( 2.0 * v1 ) * ( 2.0 - 2.0 * v1 + pow( v1, 2 ) ) ) ) / pow( Lambda[1], 3 );
    /*
        y( 0, 1 ) = ( exp( Lambda[0] - Lambda[1] ) * ( 1.0 + exp( 2.0 * Lambda[1] ) * ( -1.0 + Lambda[1] ) + Lambda[1] ) ) / pow( Lambda[1], 2 );
        y( 1, 0 ) = y( 0, 1 );
        y( 1, 1 ) = ( exp( Lambda[0] - Lambda[1] ) *
                      ( -2.0 - 2.0 * Lambda[1] - pow( Lambda[1], 2 ) + exp( 2.0 * Lambda[1] ) * ( 2.0 - 2.0 * Lambda[1] + pow( Lambda[1], 2 ) ) )
       ) / pow( Lambda[1], 3 ); */
    // std::cout << y << std::endl;
}

void M1IPMClosure::DS( Vector& ds, const Vector& u ) const {
    double alpha1 = Bisection( -10.0, 10.1, u[1] / u[0] );
    double alpha0 = -log( ( exp( alpha1 ) - exp( -alpha1 ) ) / ( alpha1 * u[0] ) );

    ds[0] = alpha0;
    ds[1] = alpha1;
    std::cout << ds << std::endl;
}

double M1IPMClosure::RootFun( const double alpha, const double u1Du0 ) const { return MathTools::coth( alpha ) - 1.0 / alpha - u1Du0; }

double M1IPMClosure::Bisection( double alphaA, double alphaB, const double u1Du0 ) const {
    if( RootFun( alphaA, u1Du0 ) * RootFun( alphaB, u1Du0 ) >= 0 ) {
        std::cout << "Incorrect a and b with vals " << RootFun( alphaA, u1Du0 ) << ", " << RootFun( alphaB, u1Du0 ) << std::endl;
        exit( EXIT_FAILURE );
    }

    double alphaC = alphaA;
    double e      = 1e-15;

    while( std::fabs( alphaB - alphaA ) >= e ) {
        alphaC = ( alphaA + alphaB ) / 2;
        if( !std::isfinite( RootFun( alphaC, u1Du0 ) ) ) {
            std::cerr << "[Bisection] Infinite RootFunction" << std::endl;
            exit( EXIT_FAILURE );
        }
        else if( RootFun( alphaC, u1Du0 ) == 0.0 ) {
            return alphaC;
        }
        else if( RootFun( alphaC, u1Du0 ) * RootFun( alphaA, u1Du0 ) < 0 ) {
            alphaB = alphaC;
        }
        else {
            alphaA = alphaC;
        }
    }
    std::cout << "result = " << alphaC << std::endl;
    // exit( EXIT_FAILURE );
    return alphaC;
}
