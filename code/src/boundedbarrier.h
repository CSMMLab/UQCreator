#ifndef BOUNDEDBARRIER_H
#define BOUNDEDBARRIER_H

#include "closure.h"
#include "problem.h"

class BoundedBarrier : public Closure
{
  private:
    BoundedBarrier() = delete;

  public:
    BoundedBarrier( Settings* settings );
    virtual ~BoundedBarrier();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif    // BOUNDEDBARRIER_H
