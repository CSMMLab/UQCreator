#include "legendre.h"

Legendre::Legendre( unsigned degree ) : Polynomial( degree ) {
    _nodes   = Vector( _degree, 0.0 );
    _weights = Vector( _degree, 0.0 );
    Compute();
}

void Legendre::Compute() {
    // construct companion matrix
    Matrix CM( _degree, _degree, 0.0 );

    for( unsigned i = 0; i < _degree - 1; ++i ) {
        CM( i + 1, i ) = std::sqrt( 1 / ( 4 - 1 / std::pow( static_cast<double>( i + 1 ), 2 ) ) );
        CM( i, i + 1 ) = std::sqrt( 1 / ( 4 - 1 / std::pow( static_cast<double>( i + 1 ), 2 ) ) );
    }

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix( CM );

    for( unsigned i = 0; i < _degree; ++i ) {
        if( std::fabs( evSys.first[i] ) < 1e-15 )
            _nodes[i] = 0;
        else
            _nodes[i] = evSys.first[i];
        _weights[i] = 2 * std::pow( evSys.second( 0, i ), 2 );
    }
    Sort();
}

double Legendre::Evaluate( unsigned m, double x ) { return boost::math::legendre_p( m, x ); }

const Vector& Legendre::GetNodes() { return _nodes; }

const Vector& Legendre::GetWeights() { return _weights; }

double Legendre::fXi( const double xi ) const { return 0.5; }

double Legendre::L2NormSquare( unsigned i ) const { return 1.0 / ( 2.0 * i + 1.0 ); }
