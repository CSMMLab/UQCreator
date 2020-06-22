#ifndef EULERCLOSURE_H
#define EULERCLOSURE_H

#include "closure.h"

class EulerClosure : public Closure
{
  private:
    double _gamma;
    EulerClosure() = delete;

  public:
    EulerClosure( Settings* settings );
    virtual ~EulerClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif    // EULERCLOSURE_H
