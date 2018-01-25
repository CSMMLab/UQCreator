#ifndef STOCHASTICGALERKIN_H
#define STOCHASTICGALERKIN_H

#include "closure.h"

class StochasticGalerkin : public Closure
{
public:
    StochasticGalerkin() = delete;
    StochasticGalerkin(Problem* problem);
    virtual ~StochasticGalerkin();

    virtual void U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda );
    virtual blaze::DynamicMatrix<double> U( const blaze::DynamicMatrix<double>& Lambda );
    virtual void DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda );
};

#endif // STOCHASTICGALERKIN_H
