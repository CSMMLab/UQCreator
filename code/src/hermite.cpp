#include "hermite.h"

Hermite::Hermite( unsigned degree ) : Polynomial( degree ) {
    _nodes.resize( degree );
    _weights.resize( degree );
    Compute();
}

void Hermite::Compute() {
    assert( _degree > 0 );

    // construct companion matrix
    Matrix CM( _degree, _degree );

    for( unsigned i = 0; i < _degree - 1; ++i ) {
        CM( i + 1, i ) = std::sqrt( static_cast<double>( i + 1 ) / 2.0 );
        CM( i, i + 1 ) = std::sqrt( static_cast<double>( i + 1 ) / 2.0 );
    }

    std::cout << CM << std::endl;

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix( CM );

    for( unsigned i = 0; i < _degree; ++i ) {
        _nodes[i]   = evSys.first[i];
        _weights[i] = std::pow( evSys.second( 0, i ), 2 ) * std::sqrt( PI );
    }
    Sort();
}

double Hermite::Evaluate( unsigned m, double x ) { return boost::math::hermite( m, x ); }

const Vector& Hermite::GetNodes() { return _nodes; }

const Vector& Hermite::GetWeights() { return _weights; }

double Hermite::fXi( const double xi ) const {
    return 1.0;    // * exp( -0.5 * pow( xi, 2 ) );
}

double Hermite::L2NormSquare( unsigned i ) const { return factorial( i ) * sqrt( PI ) * pow( 2, i ); }

unsigned Hermite::factorial( unsigned n ) const { return ( n == 1 || n == 0 ) ? 1 : factorial( n - 1 ) * n; }
