#ifndef NOLIMITER_H
#define NOLIMITER_H

#include "limiter.h"
#include "closure.h"
#include "problem.h"

class NoLimiter : public Limiter
{

public:
    NoLimiter(Closure* pClosure, Problem *problem);
    double CalculateSlope(double u0, double u1, double u2);
};

#endif // NOLIMITER_H
