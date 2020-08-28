#include "stochasticgalerkin.h"

StochasticGalerkin::StochasticGalerkin( Settings* settings ) : Closure( settings ) { _alpha = 1.0; }

StochasticGalerkin::~StochasticGalerkin() {}

void StochasticGalerkin::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void StochasticGalerkin::U( Tensor& out, const Tensor& Lambda ) { out = Lambda; }

Tensor StochasticGalerkin::U( const Tensor& Lambda ) { return Lambda; }

void StochasticGalerkin::DU( Matrix& y, const Vector& Lambda ) {
    unused( Lambda );

    y = VectorSpace::IdentityMatrix<double>( _nStates );
}

void StochasticGalerkin::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    unused( refLevel );

    _filter->FilterMoments( lambda, u );
}

void StochasticGalerkin::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    unused( refLevel );

    _filter->FilterMoments( lambda, u );
}
