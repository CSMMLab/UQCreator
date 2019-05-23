#ifndef REGULARIZEDEXPEULER1D_H
#define REGULARIZEDEXPEULER1D_H

#include "eulerclosure.h"

class RegularizedExpEuler1D : public EulerClosure
{
  private:
    double _gamma;
    double _eta;       // regularization parameter
    double _lambda;    // filter strength
    Vector _filterFunction;
    RegularizedExpEuler1D() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual void GradientNoRegularization( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );

  public:
    RegularizedExpEuler1D( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual ~RegularizedExpEuler1D();
};

#endif // REGULARIZEDEXPEULER1D_H

