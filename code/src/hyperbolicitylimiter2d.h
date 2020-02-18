#ifndef HYPERBOLICITYLIMITER2D_H
#define HYPERBOLICITYLIMITER2D_H

#include "closure.h"

class HyperbolicityLimiter2D : public Closure
{
  private:
    HyperbolicityLimiter2D() = delete;
    double _gamma;
    Vector _filterFunction;
    double _lambda;

  public:
    HyperbolicityLimiter2D( Settings* settings );
    virtual ~HyperbolicityLimiter2D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
};

#endif // HYPERBOLICITYLIMITER2D_H
