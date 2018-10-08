#include "euler.h"

Euler::Euler( std::string inputFile ) : Problem( inputFile ) {
    try {
        auto file = cpptoml::parse_file( _inputFile );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

double Euler::GetGamma() const { return _gamma; }

void Euler::Solve() {}

void Euler::Print() {}

void Euler::WriteToFile( std::string filename, int filetype ) const {}

Vector Euler::G( const Vector& u, const Vector& v ) {
    double rhoInv = 1.0 / u[1];
    double vU     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[1] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[1];
    double vV     = u[2] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[3] - 0.5 * v[1] * pow( vV, 2 ) );
    double aV     = sqrt( _gamma * p * rhoInv );

    double lambdaMin = vU - aU;
    double lambdaMax = vV + aV;

    if( lambdaMin >= 0 )
        return F( u );
    else if( lambdaMax <= 0 )
        return F( v );
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * F( u ) - lambdaMin * F( v ) + lambdaMax * lambdaMin * ( v - u ) );
    }
}

Matrix Euler::G( const Matrix& u, const Matrix& v ) {
    Matrix y( u.columns(), u.rows() );
    Vector uVec( u.columns() );
    Vector vVec( v.columns() );
    for( unsigned k = 0; k < u.rows(); ++k ) {
        for( unsigned l = 0; l < u.columns(); ++l ) {
            uVec[l] = u( l, k );
            vVec[l] = v( l, k );
        }
        Vector out = G( uVec, vVec );
        for( unsigned l = 0; l < u.columns(); ++l ) {
            y( l, k ) = out[l];
        }
    }
    return y;
}

double Euler::ExactSolution( double t, double x, double xi ) {}

Vector Euler::F( const Vector& u ) {
    double rhoInv = 1.0 / u[1];
    double v      = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[1] * pow( v, 2 ) );
    Vector flux( u.size() );
    flux[1] = u[2];
    flux[2] = u[2] * v + p;
    flux[3] = ( u[3] + p ) * v;
    return flux;
}

Matrix Euler::F( const Matrix& u ) { return 0.5 * blaze::pow( u, 2 ); }
