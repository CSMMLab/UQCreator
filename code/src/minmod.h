#ifndef MINMOD_H
#define MINMOD_H

#include "limiter.h"
#include "closure.h"
#include "problem.h"

class Minmod : public Limiter
{
private:
    virtual double CalculateSlope(double u0, double u1, double u2);
    double minmod(double a, double b);
public:
    Minmod(Closure* pClosure, Problem* problem);
};

#endif // MINMOD_H
