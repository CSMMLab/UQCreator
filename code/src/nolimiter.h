#ifndef NOLIMITER_H
#define NOLIMITER_H

#include "limiter.h"

class NoLimiter : public Limiter
{

  public:
    NoLimiter( Settings* settings, Mesh* mesh, Closure* closure );
    virtual ~NoLimiter();
    virtual double CalculateSlope( const double& u0, const double& u1, const double& u2 );
};

#endif    // NOLIMITER_H
