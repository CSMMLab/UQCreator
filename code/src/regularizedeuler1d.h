#ifndef REGULARIZEDEULER1D_H
#define REGULARIZEDEULER1D_H

#include "eulerclosure.h"

class RegularizedEuler1D : public EulerClosure
{
  private:
    double _gamma;
    double _eta;       // regularization parameter
    double _lambda;    // filter strength
    Vector _filterFunction;
    RegularizedEuler1D() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual void GradientNoRegularization( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );

  public:
    RegularizedEuler1D( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual ~RegularizedEuler1D();
};

#endif    // REGULARIZEDEULER1D_H
