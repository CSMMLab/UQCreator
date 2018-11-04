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
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif    // EULERCLOSURE_H
