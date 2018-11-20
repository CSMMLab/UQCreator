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
