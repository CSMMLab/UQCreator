#include "EulerEquationsIdGas.h"

#include <cmath>

EulerEquationsIdGas::EulerEquationsIdGas( Settings* settings ) : PhysicalProblem( settings ) {
    if( _settings->HeatingOn() ) {
        _q0       = _settings->GetQ0();
        _maxWidth = _settings->GetMaxHeatingWidth();
    }
}

EulerEquationsIdGas::~EulerEquationsIdGas() {}

void EulerEquationsIdGas::AnalyticFlux( Vector* state, Vector* flux ) {

    TransformToPrimitives( state );
    ( *flux )[0] = ( *state )[1];
    ( *flux )[1] = ( *state )[1] * _u + _p;
    ( *flux )[2] = ( ( *state )[2] + _p ) * _u;
}

void EulerEquationsIdGas::EstimateMinSpeed( Vector* state, double& minSpeed ) {
    TransformToPrimitives( state );
    minSpeed = ( _u - _a );
}

void EulerEquationsIdGas::t1_EstimateMinSpeed( Vector* state, Vector* t1_state, double& minSpeed, double& t1_minSpeed ) {
    t1_TransformToPrimitives( state, t1_state );
    TransformToPrimitives( state );

    t1_minSpeed = _t1_u + _t1_a;
    minSpeed    = _u - _a;
}

void EulerEquationsIdGas::EstimateMaxSpeed( Vector* state, double& maxSpeed ) {
    TransformToPrimitives( state );
    maxSpeed = ( _u + _a );
}

void EulerEquationsIdGas::t1_EstimateMaxSpeed( Vector* state, Vector* t1_state, double& maxSpeed, double& t1_maxSpeed ) {
    t1_TransformToPrimitives( state, t1_state );
    TransformToPrimitives( state );

    t1_maxSpeed = _t1_u + _t1_a;
    maxSpeed    = _u + _a;
}

double EulerEquationsIdGas::CFL( Vector* state, double dt, double dx ) {
    double speed, speedMax, speedMin, cfl;

    EstimateMaxSpeed( state, speedMax );
    EstimateMinSpeed( state, speedMin );
    speed = std::max( fabs( speedMax ), fabs( speedMin ) );
    cfl   = speed * dt / dx;
    if( cfl > _cfl ) _cfl = cfl;
    return speed * dt / dx;
}

void EulerEquationsIdGas::TransformToPrimitives( Vector* state ) {
    double kineticEnergy = 0;
    double innerEnergy   = 0;
    double rhoInv;
    _rho          = ( *state )[0];
    rhoInv        = 1 / _rho;
    _u            = ( *state )[1] * rhoInv;
    kineticEnergy = 0.5 * _rho * pow( _u, 2 );
    innerEnergy   = ( ( *state )[2] - kineticEnergy ) * rhoInv;
    _p            = _rho * ( _Gamma - 1 ) * innerEnergy;
    _a            = sqrt( _Gamma * _p * rhoInv );
    _T            = _p * rhoInv / _SpecificR;
}

void EulerEquationsIdGas::TransformToConservatives( double rho, double p, double u, Vector* state ) {
    double kineticEnergy = 0;
    double innerEnergy   = 0;

    ( *state )[0] = rho;
    ( *state )[1] = rho * u;
    kineticEnergy = 0.5 * rho * pow( u, 2 );
    innerEnergy   = ( p / ( rho * ( _Gamma - 1 ) ) ) * rho;
    ( *state )[2] = kineticEnergy + innerEnergy;
}
void EulerEquationsIdGas::GetPrimitives( Vector* state, double& p, double& u ) {
    TransformToPrimitives( state );
    p = _p;
    u = _u;
}
void EulerEquationsIdGas::GetMisc( Vector* state, double& T, double& mach ) {
    TransformToPrimitives( state );
    T    = _T;
    mach = std::abs( _u ) / _a;
}
double EulerEquationsIdGas::GetPressure( Vector* state ) {
    TransformToPrimitives( state );
    return _p;
}

double EulerEquationsIdGas::GetVelocity( Vector* state ) {
    TransformToPrimitives( state );
    return _u;
}
double EulerEquationsIdGas::GetTemperature( Vector* state ) {
    TransformToPrimitives( state );
    return _T;
}
double EulerEquationsIdGas::GetMachNumber( Vector* state ) {
    TransformToPrimitives( state );
    return _u / _a;
}

unsigned int EulerEquationsIdGas::GetStateDim() const { return _StateDim; }

double EulerEquationsIdGas::HeatFlux( double x ) const { return _q0 * exp( -x * x / _maxWidth ); }
