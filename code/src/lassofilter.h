#ifndef LASSOFILTER_H
#define LASSOFILTER_H

#include "closure.h"

class LassoFilter : public Closure
{
  private:
    LassoFilter() = delete;
    double _lambda;
    Vector _filterParam;
    Vector _l1Norms;

  public:
    LassoFilter( Settings* settings );
    virtual ~LassoFilter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u );
};

#endif // LASSOFILTER_H
