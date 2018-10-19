#include "burgers.h"

Burgers::Burgers( const Settings* settings ) : Problem( settings ) {
    _nStates = 1;
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        // TODO
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _settings->GetInputFile() << ": " << e.what() << std::endl;
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
    unsigned long nStates = u.rows();
    unsigned long Nq      = u.columns();
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

Matrix Burgers::F( const Matrix& u ) { return 0.5 * blaze::pow( u, 2 ); }

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

double Burgers::ExactSolution( double t, double x, double xi ) {
    double uL    = 12.0;
    double uR    = 3.0;
    double sigma = 0.2;
    double x0    = 0.5;
    double x1    = 1.5;
    double x0Bar;
    double x1Bar;
    bool shock = true;
    double y   = 0;
    if( shock ) {
        x0Bar = x0 + sigma * xi + uL * ( 1.0 / 9.0 ) + 0.5 * ( uL + uR ) * ( t - ( 1.0 / 9.0 ) );
        x1Bar = x1 + sigma * xi + uR * ( 1.0 / 9.0 ) + 0.5 * ( uL + uR ) * ( t - ( 1.0 / 9.0 ) );
        if( x <= x0Bar ) {
            return uL;
        }
        else if( x > x1Bar ) {
            return uR;
        }
    }
    else {
        x0Bar = x0 + sigma * xi + uL * t;
        x1Bar = x1 + sigma * xi + uR * t;
        if( x < x0Bar ) {
            y = uL;
        }
        else if( x > x1Bar ) {
            y = uR;
        }
        else {
            y = uL + ( uR - uL ) * ( x0Bar - x ) / ( x0Bar - x1Bar );
        }
    }
    return y;
}
