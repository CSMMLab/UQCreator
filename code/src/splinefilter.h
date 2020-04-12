#ifndef SPLINEFILTER_H
#define SPLINEFILTER_H

#include "closure.h"

class SplineFilter : public Closure
{
  private:
    SplineFilter() = delete;
    double _lambda; // filter strength
    Vector _filterFunction;
    unsigned _filterOrder;    // order of the filter
    double _eta; // variance correction

    double FilterFunction(double eta)const;

  public:
    SplineFilter( Settings* settings );
    virtual ~SplineFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
};

#endif // SPLINEFILTER_H
