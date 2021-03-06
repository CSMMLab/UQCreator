#ifndef L2FILTER_H
#define L2FILTER_H

#include "closure.h"

class L2Filter : public Closure
{
  private:
    L2Filter() = delete;
    double _lambda;
    Vector _filterFunction;

  public:
    L2Filter( Settings* settings );
    virtual ~L2Filter();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
};

#endif    // L2FILTER_H
