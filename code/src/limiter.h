#ifndef LIMITER_H
#define LIMITER_H

#include "blaze/math/DynamicVector.h"
#include "closure.h"
#include "problem.h"

class Limiter
{
private:
    Closure* _closure;
    blaze::DynamicVector<double> SlopeInternal(blaze::DynamicVector<double> u0, blaze::DynamicVector<double> u1, blaze::DynamicVector<double> u2);
    double SlopeBoundPres(double u, double slope);
protected:
    Problem* _problem;
public:
    Limiter(Closure* pClosure, Problem* problem);
    blaze::DynamicVector<double> Slope(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2, const blaze::DynamicVector<double>& lambda3);
    virtual double CalculateSlope(double u0, double u1, double u2)=0;
private:
    Limiter(){}
  //virtual ~Limiter();
};

#endif // LIMITER_H
