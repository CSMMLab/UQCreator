#include "euler.h"

Euler::Euler( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _settings->GetInputFile() << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

Euler::~Euler() {}

void Euler::Solve() {}

Vector Euler::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    double rhoInv = 1.0 / u[0];
    double vU     = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[0];
    double vV     = v[1] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[2] - 0.5 * v[0] * pow( vV, 2 ) );
    double aV     = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * vU;
    double uVProjected = nUnit[0] * vV;

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

Matrix Euler::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix Euler::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v + p;
    flux( 2, 0 ) = ( u[2] + p ) * v;
    return flux;
}

Matrix Euler::F( const Matrix& u ) {
    std::cerr << "Flux not implemented" << std::endl;
    exit( EXIT_FAILURE );
}

MatVec Euler::InitLambda( const MatVec& u ) {
    MatVec lambda( _settings->GetNumCells() + 1, Matrix( _settings->GetNStates(), _settings->GetNMoments(), 0.0 ) );
    for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
        double gamma      = -_settings->GetGamma();
        double rho        = u[j]( 0, 0 );
        double rhoU       = u[j]( 1, 0 );
        double rhoU2      = pow( rhoU, 2 );
        double rhoE       = u[j]( 2, 0 );
        lambda[j]( 0, 0 ) = ( rhoU2 + gamma * ( 2 * rho * rhoE - rhoU2 ) ) / ( -2 * rho * rhoE + rhoU2 ) -
                            std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 ) / ( 2 * rho ) ) );
        lambda[j]( 1, 0 ) = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 ) );
        lambda[j]( 2, 0 ) = -( rho / ( rhoE - ( rhoU2 ) / ( 2 * rho ) ) );
    }
    return lambda;
}

double Euler::ComputeDt( Vector& u, double dx ) const {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    double a      = sqrt( _gamma * p * rhoInv );

    double dt1 = dx * _settings->GetCFL() * ( v - a );
    double dt2 = dx * _settings->GetCFL() * ( v + a );
    if( dt1 < dt2 )
        return dt1;
    else
        return dt2;
}
