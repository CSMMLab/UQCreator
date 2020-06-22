#ifndef BOUNDEDBARRIER_H
#define BOUNDEDBARRIER_H

#include "closure.h"
#include "problem.h"

class BoundedBarrier : public Closure
{
  private:
    double _uMinus, _uPlus;
    BoundedBarrier() = delete;

  public:
    BoundedBarrier( Settings* settings );
    virtual ~BoundedBarrier();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif    // BOUNDEDBARRIER_H
