#ifndef M1IPMCLOSURE_H
#define M1IPMCLOSURE_H

#include "closure.h"

class M1IPMClosure : public Closure
{
  private:
    double _gamma;
    M1IPMClosure() = delete;

    double RootFun( const double alpha, const double u1Du0 ) const;
    double Bisection( double alphaA, double alphaB, const double u1Du0 ) const;

  public:
    M1IPMClosure( Settings* settings );
    virtual ~M1IPMClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif    // M1IPMCLOSURE_H
