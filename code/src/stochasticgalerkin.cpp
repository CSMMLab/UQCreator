#include "stochasticgalerkin.h"

StochasticGalerkin::StochasticGalerkin( Settings* settings ) : Closure( settings ) { _alpha = 1.0; }

StochasticGalerkin::~StochasticGalerkin() {}

void StochasticGalerkin::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void StochasticGalerkin::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix StochasticGalerkin::U( const Matrix& Lambda ) { return Lambda; }

void StochasticGalerkin::DU( Matrix& y, const Vector& Lambda ) { y = blaze::IdentityMatrix<double>( _nStates ); }

void StochasticGalerkin::SolveClosure( Matrix& lambda, const Matrix& u ) { lambda = u; }
