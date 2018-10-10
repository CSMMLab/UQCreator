#include "Advection.h"

#include <cmath>

Advection::Advection( Settings* settings ) : PhysicalProblem( settings ) {}

Advection::~Advection() {}

void Advection::AnalyticFlux( Vector* state, Vector* flux ) { ( *flux )[0] = _a * ( *state )[0]; }

void Advection::EstimateMinSpeed( Vector* state, double& minSpeed ) { minSpeed = _a; }

void Advection::EstimateMaxSpeed( Vector* state, double& maxSpeed ) { maxSpeed = _a; }

unsigned int Advection::GetStateDim() const { return _StateDim; }

double Advection::CFL( Vector* state, double dt, double dx ) { return _a * dt / dx; }
