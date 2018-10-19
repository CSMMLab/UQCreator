#ifndef NOLIMITER_H
#define NOLIMITER_H

#include "closure.h"
#include "limiter.h"
#include "problem.h"

class NoLimiter : public Limiter
{

  public:
    NoLimiter( const Closure* closure, const Settings* settings );
    virtual ~NoLimiter();
    virtual double CalculateSlope( const double& u0, const double& u1, const double& u2 );
};

#endif    // NOLIMITER_H
