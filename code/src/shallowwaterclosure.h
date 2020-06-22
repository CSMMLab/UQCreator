#ifndef SHALLOWWATERCLOSURE_H
#define SHALLOWWATERCLOSURE_H

#include "closure.h"
#include "problem.h"

class ShallowWaterClosure : public Closure
{
  private:
    ShallowWaterClosure() = delete;
    double _g;

  public:
    ShallowWaterClosure( Settings* settings );
    virtual ~ShallowWaterClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif    // SHALLOWWATERCLOSURE_H
