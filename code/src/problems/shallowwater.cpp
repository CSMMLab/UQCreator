#include "problems/shallowwater.h"

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

Matrix ShallowWater::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
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

Matrix ShallowWater::F( const Matrix& /*u*/ ) {
    _log->error( "[ShallowWater] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ShallowWater::ComputeDt( const Matrix& /*u*/, double /*dx*/, unsigned /*level*/ ) const {
    _log->error( "[ShallowWater] ComputeDt not implemented" );
    return 0.0;
}

Vector ShallowWater::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );
    _sigma = _settings->GetSigma();
    if( xi.size() == 1 ) {
        double a    = 1000.0;
        double uL   = 10.0;
        double hLhR = 0.05;
        double uR   = uL * hLhR;
        y[1]        = 0.0;
        if( x[0] < a + _sigma[0] * xi[0] ) {
            y[0] = uL;
            return y;
        }
        else {
            y[0] = uR;
            return y;
        }
    }
    else if( xi.size() == 2 ) {
        double x0 = 0.3;
        double x1 = 0.6;
        double uL = 12.0;
        double uM = 6.0;
        double uR = 1.0;

        if( x[0] < x0 )
            y[0] = uL + _sigma[0] * xi[0];
        else if( x[0] < x1 )
            y[0] = uM + _sigma[1] * xi[1];
        else
            y[0] = uR;
        return y;
    }
    _log->error( "Reached end of IC. No initial condition set" );
    exit( EXIT_FAILURE );
}

Vector ShallowWater::LoadIC( const Vector& /*x*/, const Vector& /*xi*/ ) {
    _log->error( "[ShallowWater: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
