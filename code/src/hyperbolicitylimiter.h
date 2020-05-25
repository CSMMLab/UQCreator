#ifndef HYPERBOLICITYLIMITER_H
#define HYPERBOLICITYLIMITER_H

#include "closure.h"

class HyperbolicityLimiter : public Closure
{
  private:
    HyperbolicityLimiter() = delete;
    double _gamma;
    Vector _filterFunction;
    unsigned _filterOrder;    // order of the filter
    double _lambda;
    double _c;    // machine precision constant
    double FilterFunction( double eta ) const;

  public:
    HyperbolicityLimiter( Settings* settings );
    virtual ~HyperbolicityLimiter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
};

#endif    // HYPERBOLICITYLIMITER_H
