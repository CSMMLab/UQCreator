#include "advection2d.h"

Advection2D::Advection2D( Settings* settings ) : Problem( settings ) {
    _nStates = 1;
    _settings->SetNStates( _nStates );
    _omega = Vector( 2, 1.0 ) / sqrt( 2.0 );

    // for the one dimensional case an exact solution is specified
    _settings->SetExactSolution( false );

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Advection2D] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Vector Advection2D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    unused( nUnit );

    double inner = _omega[0] * n[0] + _omega[1] * n[1];
    if( inner > 0 ) {
        return inner * u;
    }
    else {
        return inner * v;
    }
}

Matrix Advection2D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Vector Advection2D::F( double u ) {
    Vector y( 2 );
    y[0] = _omega[0] * u;
    y[1] = _omega[1] * u;
    return y;
}

Vector Advection2D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );
    _sigma = _settings->GetSigma();

    Vector xM( 2, 0.75 + _sigma[0] * xi[0] );

    double t = 0.01;    // pseudo time for gaussian smoothing
    y[0]     = 1.0 / ( 4.0 * M_PI * t ) * std::exp( -( ( x[0] - xM[0] ) * ( x[0] - xM[0] ) + ( x[1] - xM[1] ) * ( x[1] - xM[1] ) ) / ( 4 * t ) );
    return y;
}

Matrix Advection2D::ExactSolution( double t, const Matrix& x, const Vector& xi ) const {
    double x0, x1;

    Matrix y( _settings->GetNumCells(), _nStates );
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

    for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
        if( x( j, 0 ) < x0 )
            y( j, 0 ) = _uL;
        else if( x( j, 0 ) < x1 )
            y( j, 0 ) = _uL + ( _uR - _uL ) * ( x( j, 0 ) - x0 ) / ( x1 - x0 );
        else
            y( j, 0 ) = _uR;
    }
    return y;
}

double Advection2D::ComputeDt( const Tensor& u, double dx, unsigned level ) const {
    unused( u );
    unused( level );

    return _settings->GetCFL() / dx / norm( _omega );
}

Vector Advection2D::LoadIC( const Vector& x, const Vector& xi ) {
    unused( x );
    unused( xi );

    _log->error( "[Advection2D: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}
