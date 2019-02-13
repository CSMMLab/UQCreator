#include "burgers.h"

Burgers::Burgers( Settings* settings ) : Problem( settings ) {
    _nStates = 1;
    _settings->SetNStates( _nStates );

    // for the one dimensional case an exact solution is specified
    if( _settings->GetNDimXi() == 1 ) {
        _settings->SetExactSolution( true );
    }

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[burgers] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
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
    unsigned nStates = u.rows();
    unsigned Nq      = u.columns();
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

Matrix Burgers::F( const Matrix& u ) { return 0.5 * pow( u, 2 ); }

Vector Burgers::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );
    _sigma = Vector( xi.size() );
    if( xi.size() == 1 ) {

        _sigma[0] = 0.1;    // 0.2
        _x0       = 0.5;
        _x1       = 1.5;
        _uL       = 12.0;
        _uR       = 3.0;
        if( x[0] < _x0 + _sigma[0] * xi[0] ) {
            y[0] = _uL;
            return y;
        }
        else if( x[0] < _x1 + _sigma[0] * xi[0] ) {
            y[0] = _uL + ( _uR - _uL ) * ( _x0 + _sigma[0] * xi[0] - x[0] ) / ( _x0 - _x1 );
            return y;
        }
        else {
            y[0] = _uR;
            return y;
        }
    }
    else if( xi.size() == 2 ) {
        _x0           = 0.3;
        _x1           = 0.6;
        double sigma0 = 0.2;    // 0.2
        double sigma1 = 0.1;
        _uL           = 12.0;
        double uM     = 6.0;
        _uR           = 1.0;

        if( x[0] < _x0 )
            y[0] = _uL + sigma0 * xi[0];
        else if( x[0] < _x1 )
            y[0] = uM + sigma1 * xi[1];
        else
            y[0] = _uR;
        return y;
    }
    _log->error( "Reached end of IC. No initial condition set" );
    exit( EXIT_FAILURE );
}

Vector Burgers::ExactSolution( double t, const Vector& x, const Vector& xi ) const {
    double x0, x1;
    Vector y( _nStates );
    if( t >= ( _x1 - _x0 ) / ( _uL - _uR ) ) {
        double tS            = ( _x1 - _x0 ) / ( _uL - _uR );
        double x0BeforeShock = _x0 + _sigma[0] * xi[0] + tS * _uL;
        // double x1BeforeShock = _x1 + _sigma[0] * xi[0] + tS * _uR;
        x0 = x0BeforeShock + ( t - tS ) * ( _uL + _uR ) * 0.5;
        x1 = x0 - 1.0;
    }
    else {
        x0 = _x0 + _sigma[0] * xi[0] + t * _uL;
        x1 = _x1 + _sigma[0] * xi[0] + t * _uR;
    }

    if( x[0] < x0 )
        y[0] = _uL;
    else if( x[0] < x1 )
        y[0] = _uL + ( _uR - _uL ) * ( x[0] - x0 ) / ( x1 - x0 );
    else
        y[0] = _uR;
    return y;
}

double Burgers::ComputeDt( const Matrix& u, double dx ) const { return dx * _settings->GetCFL() / u( 0, 0 ); }

Vector Burgers::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[Burgers: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
