#ifndef REGULARIZEDLASSOEULER_H
#define REGULARIZEDLASSOEULER_H

#include "eulerclosure.h"

class RegularizedLassoEuler : public EulerClosure
{
  private:
    double _gamma;
    double _eta;       // regularization parameter
    double _lambda;    // filter strength

    Vector _filterParam;
    Vector _l1Norms;

    Vector _filterFunction;
    RegularizedLassoEuler() = delete;
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );

  public:
    RegularizedLassoEuler( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
    virtual ~RegularizedLassoEuler();
};

#endif    // REGULARIZEDLASSOEULER_H
