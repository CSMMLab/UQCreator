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

    virtual void U( Vector& out, const Vector& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // BOUNDEDBARRIER_H
