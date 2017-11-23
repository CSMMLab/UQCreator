#ifndef LIMITER_H
#define LIMITER_H

#include "blaze/math/DynamicVector.h"
#include <cpptoml.h>
#include "closure.h"
#include "problem.h"

class Limiter
{
private:
    Closure* _closure;
    blaze::DynamicVector<double> SlopeInternal(const blaze::DynamicVector<double>& u0, const blaze::DynamicVector<double>& u1, const blaze::DynamicVector<double>& u2);
    double SlopeBoundPres(const double& u, const double& slope);
protected:
    Problem* _problem;
    double _dx;
public:
    Limiter(Closure* pClosure, Problem* problem);
    virtual ~Limiter();
    static Limiter* Create(Closure* closure, Problem* problem);
    blaze::DynamicVector<double> Slope(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2, const blaze::DynamicVector<double>& lambda3);
    virtual double CalculateSlope(const double& u0, const double& u1, const double& u2)=0;
private:
    Limiter(){}
};

#endif // LIMITER_H
