#ifndef SHALLOWWATERCLOSURE2D_H
#define SHALLOWWATERCLOSURE2D_H

#include "closure.h"
#include "problem.h"

class ShallowWaterClosure2D : public Closure
{
  private:
    ShallowWaterClosure2D() = delete;
    double _g;

  public:
    ShallowWaterClosure2D( Settings* settings );
    virtual ~ShallowWaterClosure2D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif    // SHALLOWWATERCLOSURE2D_H
