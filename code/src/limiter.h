#ifndef LIMITER_H
#define LIMITER_H

#include "blaze/math/DynamicVector.h"
#include "closure.h"

class Limiter
{
private:
    Closure* _closure;
    blaze::DynamicVector<double> SlopeInternal(blaze::DynamicVector<double> u0, blaze::DynamicVector<double> u1, blaze::DynamicVector<double> u2);
    virtual blaze::DynamicVector<double> CalculateSlope(blaze::DynamicVector<double> u0, blaze::DynamicVector<double> u1, blaze::DynamicVector<double> u2)=0;
public:
    Limiter(Closure* pClosure);
    blaze::DynamicVector<double> Limiter::Slope(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2, const blaze::DynamicVector<double>& lambda3);
private:
    Limiter(){}
    virtual ~Limiter(){}
};

#endif // LIMITER_H
