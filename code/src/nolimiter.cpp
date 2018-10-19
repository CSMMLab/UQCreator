#include "nolimiter.h"

NoLimiter::NoLimiter( Settings* settings, Mesh* mesh, Closure* closure ) : Limiter( settings, mesh, closure ) {}

NoLimiter::~NoLimiter() {}

double NoLimiter::CalculateSlope( const double& u0, const double& u1, const double& u2 ) { return 0.0; }
