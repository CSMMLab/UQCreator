#include "euler2d.h"

Euler2D::Euler2D( std::string inputFile ) : Problem( inputFile ) {
    _nStates = 4;
    try {
        auto file = cpptoml::parse_file( _inputFile );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

double Euler2D::GetGamma() const { return _gamma; }

void Euler2D::Solve() {}

Vector Euler2D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    double rhoInv = 1.0 / u[0];
    double uU     = u[1] * rhoInv;
    double vU     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
    double aU     = sqrt( _gamma * p * rhoInv );

    rhoInv    = 1.0 / v[0];
    double uV = v[1] * rhoInv;
    double vV = v[2] * rhoInv;
    p         = ( _gamma - 1.0 ) * ( v[3] - 0.5 * v[0] * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
    double aV = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
    double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        return F( u ) * n;
    else if( lambdaMax <= 0 )
        return F( v ) * n;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
}

Matrix Euler2D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = u.rows();
    unsigned Nq      = u.columns();
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

double Euler2D::ExactSolution( double t, double x, double xi ) { return 0.0; }

Matrix Euler2D::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v1     = u[1] * rhoInv;
    double v2     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( v1, 2 ) + pow( v2, 2 ) ) );
    Matrix flux( u.size(), 2 );

    flux( 1, 1 ) = u[1];
    flux( 2, 1 ) = u[1] * v1 + p;
    flux( 3, 1 ) = u[1] * v2;
    flux( 4, 1 ) = ( u[3] + p ) * v1;
    flux( 1, 2 ) = u[2];
    flux( 2, 2 ) = u[2] * v1;
    flux( 3, 2 ) = u[2] * v2 + p;
    flux( 4, 2 ) = ( u[3] + p ) * v2;

    return flux;
}

Matrix Euler2D::F( const Matrix& u ) {
    std::cerr << "Flux not implemented" << std::endl;
    exit( EXIT_FAILURE );
    return 0.5 * blaze::pow( u, 2 );
}
