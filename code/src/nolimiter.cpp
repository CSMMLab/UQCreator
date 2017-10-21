#include "nolimiter.h"

NoLimiter::NoLimiter(Closure *pClosure, Problem* problem) : Limiter(pClosure,problem)
{

}

double NoLimiter::CalculateSlope(double u0, double u1, double u2){
    return 0.0;
}
