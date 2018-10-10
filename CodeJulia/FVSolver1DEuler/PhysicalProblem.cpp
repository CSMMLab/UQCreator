#include "PhysicalProblem.h"
#include <cmath>

PhysicalProblem::PhysicalProblem( Settings* settings ) : _settings( settings ) { _cfl = 0; }
PhysicalProblem::~PhysicalProblem() {}

void PhysicalProblem::TransformToPrimitives( Vector* state ) {}

void PhysicalProblem::t1_TransformToPrimitives( Vector* state, Vector* t1_state ) {}

void PhysicalProblem::TransformToConservatives( double rho, double p, double u, Vector* state ) {}
void PhysicalProblem::GetPrimitives( Vector* state, double& p, double& u ) {}
void PhysicalProblem::GetMisc( Vector* state, double& T, double& mach, double& vapour ) {}
void PhysicalProblem::GetMisc( Vector* state, double& T, double& mach ) {}
double PhysicalProblem::GetPressure( Vector* state ) { return -1.0; }
double PhysicalProblem::GetVelocity( Vector* state ) { return -1.0; }
double PhysicalProblem::GetVapourFraction( Vector* state ) { return -1.0; }
double PhysicalProblem::GetTemperature( Vector* state ) { return -1.0; }
double PhysicalProblem::GetMachNumber( Vector* state ) { return -1.0; }
double PhysicalProblem::HeatFlux( double x ) const { return -1.0; }
