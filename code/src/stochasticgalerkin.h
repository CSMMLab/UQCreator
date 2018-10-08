#ifndef STOCHASTICGALERKIN_H
#define STOCHASTICGALERKIN_H

#include "closure.h"

class StochasticGalerkin : public Closure
{
public:
    StochasticGalerkin() = delete;
    StochasticGalerkin(Problem* problem);
    virtual ~StochasticGalerkin();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // STOCHASTICGALERKIN_H
