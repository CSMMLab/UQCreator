#ifndef HYPERBOLICITYLIMITER2D_H
#define HYPERBOLICITYLIMITER2D_H

#include "closure.h"

class HyperbolicityLimiter2D : public Closure
{
  private:
    HyperbolicityLimiter2D() = delete;
    double _gamma;

  public:
    HyperbolicityLimiter2D( Settings* settings );
    virtual ~HyperbolicityLimiter2D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual void SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel );
};

#endif    // HYPERBOLICITYLIMITER2D_H
