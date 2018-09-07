#include "stochasticgalerkin.h"

StochasticGalerkin::StochasticGalerkin( Problem* problem ) : Closure( problem ) {}

StochasticGalerkin::~StochasticGalerkin() {}

void StochasticGalerkin::U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda ) { out = Lambda; }

blaze::DynamicMatrix<double> StochasticGalerkin::U( const blaze::DynamicMatrix<double>& Lambda ) { return Lambda; }

void StochasticGalerkin::DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda ) {
    y = blaze::IdentityMatrix<double>( _nStates );
}
