#ifndef MINMOD_H
#define MINMOD_H

#include "limiter.h"
#include "closure.h"
#include "problem.h"

class Minmod : public Limiter
{
private:
    virtual double CalculateSlope(const double& u0, const double& u1, const double& u2);
    double minmod(const double& a, const double& b);
public:
    Minmod(Closure* pClosure, Problem* problem);
    virtual ~Minmod();
};

#endif // MINMOD_H