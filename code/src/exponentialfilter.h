#ifndef EXPONENTIALFILTER_H
#define EXPONENTIALFILTER_H


#include "closure.h"

class ExponentialFilter : public Closure
{
  private:
    ExponentialFilter() = delete;
    double _lambda;
    Vector _filterFunction;

  public:
    ExponentialFilter( Settings* settings );
    virtual ~ExponentialFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
};
#endif // EXPONENTIALFILTER_H
