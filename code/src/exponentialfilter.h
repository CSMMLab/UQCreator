#ifndef EXPONENTIALFILTER_H
#define EXPONENTIALFILTER_H

#include "closure.h"

class ExponentialFilter : public Closure
{
  private:
    ExponentialFilter() = delete;
    double _lambda;    // filter strength
    Vector _filterFunction;
    double _c;                // machine precision constant
    unsigned _filterOrder;    // order of the filter

    double FilterFunction( double eta ) const;

  public:
    ExponentialFilter( Settings* settings );
    virtual ~ExponentialFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Tensor );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
};
#endif    // EXPONENTIALFILTER_H
