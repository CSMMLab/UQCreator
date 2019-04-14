#ifndef REGULARIZEDBOUNDEDBARRIER_H
#define REGULARIZEDBOUNDEDBARRIER_H

#include "boundedbarrier.h"

class RegularizedBoundedBarrier : public BoundedBarrier
{
    double _eta; // regularization parameter
    double _lambda; // filter strength
    Vector _filterFunction;
    RegularizedBoundedBarrier() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
public:
    RegularizedBoundedBarrier( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual ~RegularizedBoundedBarrier();
};

#endif // REGULARIZEDBOUNDEDBARRIER_H
