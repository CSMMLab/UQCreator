#ifndef EULERCLOSURE2D_H
#define EULERCLOSURE2D_H

#include "closure.h"
#include "problem.h"

class EulerClosure2D : public Closure
{
  private:
    double _gamma;
    EulerClosure2D() = delete;

  public:
    EulerClosure2D( Settings* settings );
    virtual ~EulerClosure2D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );

    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif    // EULERCLOSURE2D_H
