#ifndef REGULARIZEDEULER_H
#define REGULARIZEDEULER_H

#include "eulerclosure2d.h"

class RegularizedEuler : public EulerClosure2D
{
  private:
    double _gamma;
    double _eta;       // regularization parameter
    double _lambda;    // filter strength
    Vector _filterFunction;
    RegularizedEuler() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual void GradientNoRegularizaton( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );

  public:
    RegularizedEuler( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual ~RegularizedEuler();
};

#endif    // REGULARIZEDEULER_H
