#ifndef STOCHASTICGALERKIN_H
#define STOCHASTICGALERKIN_H

#include "closure.h"

class StochasticGalerkin : public Closure
{
  private:
    StochasticGalerkin() = delete;

  public:
    StochasticGalerkin( Settings* settings );
    virtual ~StochasticGalerkin();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u );
};

#endif    // STOCHASTICGALERKIN_H
