#include "nolimiter.h"

NoLimiter::NoLimiter( const Closure* closure, const Settings* settings ) : Limiter( closure, problem ) {}

NoLimiter::~NoLimiter() {}

double NoLimiter::CalculateSlope( const double& u0, const double& u1, const double& u2 ) { return 0.0; }
