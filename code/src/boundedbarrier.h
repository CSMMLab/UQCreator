#ifndef BOUNDEDBARRIER_H
#define BOUNDEDBARRIER_H

#include "closure.h"
#include "problem.h"

class BoundedBarrier : public Closure
{
public:
    BoundedBarrier() = delete;
    BoundedBarrier(Problem* problem);
    virtual ~BoundedBarrier();

    virtual void U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda );
    virtual blaze::DynamicMatrix<double> U( const blaze::DynamicMatrix<double>& Lambda );
    virtual void DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda );
};

#endif // BOUNDEDBARRIER_H
