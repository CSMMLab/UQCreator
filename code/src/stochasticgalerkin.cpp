#include "stochasticgalerkin.h"

StochasticGalerkin::StochasticGalerkin( Problem* problem ) : Closure( problem ) {}

StochasticGalerkin::~StochasticGalerkin() {}

void StochasticGalerkin::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

Matrix StochasticGalerkin::U( const Matrix& Lambda ) { return Lambda; }

void StochasticGalerkin::DU( Matrix& y, const Vector& Lambda ) {
    y = blaze::IdentityMatrix<double>( _nStates );
}
