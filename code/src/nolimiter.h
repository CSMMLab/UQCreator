#ifndef NOLIMITER_H
#define NOLIMITER_H

#include "limiter.h"
#include "closure.h"
#include "problem.h"

class NoLimiter : public Limiter
{

public:
    NoLimiter(Closure* pClosure, Problem *problem);
    virtual ~NoLimiter();
    virtual double CalculateSlope(const double& u0, const double& u1, const double& u2);
};

#endif // NOLIMITER_H
