#ifndef HOULIFILTER_H
#define HOULIFILTER_H

#include "closure.h"

class HouLiFilter : public Closure
{
  private:
    HouLiFilter() = delete;
    double _lambda; // filter strength
    Vector _filterFunction;
    unsigned _filterOrder;    // order of the filter
    unsigned _gamma;

    double FilterFunction(double eta)const;

  public:
    HouLiFilter( Settings* settings );
    virtual ~HouLiFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
};

#endif // HOULIFILTER_H
