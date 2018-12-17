#include "shallowwater.h"

ShallowWater::ShallowWater( Settings* settings ) : Problem( settings ) {
    _nStates = 2;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _g           = 9.81;
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ShallowWater] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ShallowWater::~ShallowWater() {}

void ShallowWater::Solve() {}

Vector ShallowWater::G( const Vector& q_l, const Vector& q_r, const Vector& nUnit, const Vector& n ) {
    double u_r = nUnit[0] * q_r[1] / q_r[0];
    double c_r = sqrt( _g * q_r[0] );
    double u_l = nUnit[0] * q_l[1] / q_l[0];
    double c_l = sqrt( _g * q_l[0] );

    // Lax - Fridrichs
    double dtdx = _settings->GetCFL() / 12.0;
    return 0.5 * ( F( q_l ) * n + F( q_r ) * n ) - ( 0.5 / dtdx ) * ( q_r - q_l );

    // return F( q_hat ) * n;

    double lambda_l = u_l - c_l;
    double lambda_r = u_r + c_r;

    if( lambda_l >= 0 ) {
        return F( q_l ) * n;
    }
    else if( lambda_r <= 0 ) {
        return F( q_r ) * n;
    }
    else {
        return ( 1.0 / ( lambda_r - lambda_l ) ) * ( lambda_r * F( q_l ) * n - lambda_l * F( q_r ) * n + lambda_r * lambda_l * ( q_r - q_l ) );
    }
}

Matrix ShallowWater::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ShallowWater::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = pow( u[1], 2 ) / u[0] + 0.5 * _g * pow( u[0], 2 );
    return flux;
}

Matrix ShallowWater::F( const Matrix& u ) {
    _log->error( "[ShallowWater] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ShallowWater::ComputeDt( Vector& u, double dx ) const {
    _log->error( "[ShallowWater] ComputeDt not implemented" );
    return 0.0;
}
