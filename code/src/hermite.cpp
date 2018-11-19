#include "hermite.h"

Hermite::Hermite( unsigned degree ) : Polynomial( degree ) {
    _nodes.resize( degree );
    _weights.resize( degree );
    Compute();
}

void Hermite::Compute() {
    assert( _degree > 0 );

    // construct companion matrix
    Matrix CM( _degree, _degree, 0.0 );

    for( unsigned i = 0; i < _degree - 1; ++i ) {
        CM( i + 1, i ) = std::sqrt( ( i + 1 ) / 2 );
        CM( i, i + 1 ) = std::sqrt( ( i + 1 ) / 2 );
    }

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
    double pi = 3.14159265359;
    return ( 1.0 / sqrt( 2.0 * pi ) ) * exp( -0.5 * pow( xi, 2 ) );
}
