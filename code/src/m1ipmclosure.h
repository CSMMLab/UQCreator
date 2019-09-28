#ifndef M1IPMCLOSURE_H
#define M1IPMCLOSURE_H

#include "closure.h"

class M1IPMClosure : public Closure
{
  private:
    double _gamma;
    M1IPMClosure() = delete;

  public:
    M1IPMClosure( Settings* settings );
    virtual ~M1IPMClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif // M1IPMCLOSURE_H
