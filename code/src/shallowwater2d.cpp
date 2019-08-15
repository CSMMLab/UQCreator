#include "shallowwater2d.h"

ShallowWater2D::ShallowWater2D( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _g           = 9.81;
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ShallowWater2D] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ShallowWater2D::~ShallowWater2D() {}

void ShallowWater2D::Solve() {}

Vector ShallowWater2D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    double uU = u[1] / u[0];
    double vU = u[2] / u[0];
    double uV = v[1] / v[0];
    double vV = v[2] / v[0];
    double cU = sqrt( _g * u[0] );
    double cV = sqrt( _g * v[0] );

    double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
    double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

    double lambdaMin = uUProjected - cU;
    double lambdaMax = uVProjected + cV;

    if( lambdaMin >= 0 )
        return F( u ) * n;
    else if( lambdaMax <= 0 )
        return F( v ) * n;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
}

Matrix ShallowWater2D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ShallowWater2D::F( const Vector& u ) const {
    Matrix flux( u.size(), 2 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = pow( u[1], 2 ) / u[0] + 0.5 * _g * pow( u[0], 2 );
    flux( 2, 0 ) = u[1] * u[2] / u[0];

    flux( 0, 1 ) = u[2];
    flux( 1, 1 ) = u[1] * u[2] / u[0];
    flux( 2, 1 ) = pow( u[2], 2 ) / u[0] + 0.5 * _g * pow( u[0], 2 );
    return flux;
}

Matrix ShallowWater2D::F( const Matrix& u ) {
    _log->error( "[ShallowWater2D] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ShallowWater2D::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    _log->error( "[ShallowWater2D] ComputeDt not implemented" );
    return 0.0;
}

Vector ShallowWater2D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );
    _sigma    = _settings->GetSigma();
    double a  = 0.35;
    double uL = 10.0;    // 10.0;
    double uR = 5.0;
    y[1]      = 0.0;
    y[2]      = 0.0;
    if( xi.size() == 1 ) {
        if( x[0] < a + _sigma[0] * xi[0] ) {
            y[0] = uL;    // + sigma * xi[0];
            return y;
        }
        else {
            y[0] = uR;    // - sigma * xi[0];
            return y;
        }
    }
    else if( xi.size() == 2 ) {
        if( x[0] < a ) {
            y[0] = uL + _sigma[0] * xi[0];
            return y;
        }
        else {
            y[0] = uR + _sigma[1] * xi[1];
            return y;
        }
    }
    _log->error( "Reached end of IC. No initial condition set" );
    exit( EXIT_FAILURE );
}

Vector ShallowWater2D::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ShallowWater2D: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}

Matrix ShallowWater2D::BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    Vector uB( nStates );
    for( unsigned k = 0; k < Nq; ++k ) {
        Vector v( 2, 0.0 );
        v.reset();

        v[0]      = u( 1, k ) / u( 0, k );
        v[1]      = u( 2, k ) / u( 0, k );
        double vn = dot( nUnit, v );
        Vector Vn = vn * nUnit;
        Vector Vb = -Vn + v;
        // double velMagB = Vb[0] * Vb[0] + Vb[1] * Vb[1];
        // double velMag  = v[0] * v[0] + v[1] * v[1];
        double rho     = u( 0, k );
        uB[0]          = rho;
        uB[1]          = rho * ( Vb[0] );
        uB[2]          = rho * ( Vb[1] );
        column( y, k ) = F( uB ) * n;
    }
    return y;
}
