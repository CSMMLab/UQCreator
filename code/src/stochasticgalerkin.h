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
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual Matrix SolveClosure( const Matrix& u, Matrix& lambda );
};

#endif    // STOCHASTICGALERKIN_H
