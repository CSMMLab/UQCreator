#ifndef SPLINEFILTER_H
#define SPLINEFILTER_H

#include "closure.h"

class SplineFilter : public Closure
{
  private:
    SplineFilter() = delete;
    double _lambda;    // filter strength
    Vector _filterFunction;
    unsigned _filterOrder;    // order of the filter
    double _eta;              // variance correction

    double FilterFunction( double eta ) const;

  public:
    SplineFilter( Settings* settings );
    virtual ~SplineFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual void SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel );
};

#endif    // SPLINEFILTER_H
