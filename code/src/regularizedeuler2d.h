#ifndef REGULARIZEDEULER_H
#define REGULARIZEDEULER_H

#include "eulerclosure2d.h"

class RegularizedEuler2D : public EulerClosure2D
{
  private:
    double _eta;    // regularization parameter
    RegularizedEuler2D() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void GradientNoRegularizaton( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel );

  public:
    RegularizedEuler2D( Settings* settings );
    virtual void SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual void SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel );
    virtual ~RegularizedEuler2D();
};

#endif    // REGULARIZEDEULER_H
