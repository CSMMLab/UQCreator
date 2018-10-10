#include "BurgersEquation.h"

#include <cmath>

BurgersEquation::BurgersEquation( Settings* settings ) : PhysicalProblem( settings ) {}

BurgersEquation::~BurgersEquation() {}

void BurgersEquation::AnalyticFlux( Vector* state, Vector* flux ) { ( *flux )[0] = 0.5 * pow( ( *state )[0], 2 ); }

void BurgersEquation::EstimateMinSpeed( Vector* state, double& minSpeed ) { minSpeed = std::min( -0.1, ( *state )[0] ); }

void BurgersEquation::EstimateMaxSpeed( Vector* state, double& maxSpeed ) { maxSpeed = std::max( 0.1, ( *state )[0] ); }

unsigned int BurgersEquation::GetStateDim() const { return _StateDim; }

double BurgersEquation::CFL( Vector* state, double dt, double dx ) { return ( *state )[0] * dt / dx; }
