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
    double rhoInv = 1.0 / u[0];
    double vU     = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[0];
    double vV     = u[1] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[2] - 0.5 * v[0] * pow( vV, 2 ) );
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
    unsigned nStates = u.rows();
    unsigned Nq      = u.columns();
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ) );
    }
    return y;
}

double Euler::ExactSolution( double t, double x, double xi ) { return 0.0; }

Vector Euler::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Vector flux( u.size() );
    flux[0] = u[1];
    flux[1] = u[1] * v + p;
    flux[2] = ( u[2] + p ) * v;
    return flux;
}

Matrix Euler::F( const Matrix& u ) {
    std::cerr << "Flux not implemented" << std::endl;
    exit( EXIT_FAILURE );
    return 0.5 * blaze::pow( u, 2 );
}
