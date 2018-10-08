#ifndef LIMITER_H
#define LIMITER_H

#include "blaze/math/DynamicVector.h"
#include "closure.h"
#include "problem.h"
#include <cpptoml.h>

class Limiter
{
  private:
    Closure* _closure;
    Matrix
    SlopeInternal( const Matrix& u0, const Matrix& u1, const Matrix& u2 );
    double SlopeBoundPres( const double& u, const double& slope );

  protected:
    Problem* _problem;
    double _dx;

  public:
    Limiter( Closure* pClosure, Problem* problem );
    virtual ~Limiter();
    static Limiter* Create( Closure* closure, Problem* problem );
    Matrix
    Slope( const Matrix& lambda1, const Matrix& lambda2, const Matrix& lambda3 );
    virtual double CalculateSlope( const double& u0, const double& u1, const double& u2 ) = 0;

  private:
    Limiter() {}
};

#endif    // LIMITER_H
