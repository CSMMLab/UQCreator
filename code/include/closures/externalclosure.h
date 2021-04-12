#ifndef EXTERNALCLOSURE_H
#define EXTERNALCLOSURE_H

#include "closure.h"

class ExternalClosure : public Closure
{
  private:
    void ( *_U )( double*, double* );
    void ( *_DU )( double*, double* );
    void ( *_SolveClosure )( double*, double*, unsigned );

    ExternalClosure() = delete;

  public:
    ExternalClosure( void ( *U )( double*, double* ),
                     void ( *DU )( double*, double* ),
                     void ( *SolveClosure )( double*, double*, unsigned ),
                     Settings* settings );
    virtual ~ExternalClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Tensor& out, const Tensor& Lambda );
    virtual Tensor U( const Tensor& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual void SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel );
};

#endif    // EXTERNALCLOSURE_H
