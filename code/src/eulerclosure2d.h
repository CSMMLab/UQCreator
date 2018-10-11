#ifndef EULERCLOSURE2D_H
#define EULERCLOSURE2D_H

#include "closure.h"
#include "problem.h"

class EulerClosure2D : public Closure
{
  private:
    double _gamma;

  public:
    EulerClosure2D() = delete;
    EulerClosure2D( Problem* problem );
    virtual ~EulerClosure2D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif    // EULERCLOSURE2D_H
