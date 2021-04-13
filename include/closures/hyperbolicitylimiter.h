#ifndef HYPERBOLICITYLIMITER_H
#define HYPERBOLICITYLIMITER_H

#include "closure.h"

class HyperbolicityLimiter : public Closure
{
  private:
    HyperbolicityLimiter() = delete;
    double _gamma;

  public:
    HyperbolicityLimiter( Settings* settings );
    virtual ~HyperbolicityLimiter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual void SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel );
};

#endif    // HYPERBOLICITYLIMITER_H