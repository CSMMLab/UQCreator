#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "mathtools.h"
#include "polynomial.h"

#include <boost/math/special_functions/legendre.hpp>

class Legendre : public Polynomial
{
  private:
    virtual void Compute();

  public:
    Legendre( int degree );

    virtual double Evaluate( int m, double x );
    virtual const blaze::DynamicVector<double>& GetNodes();
    virtual const blaze::DynamicVector<double>& GetWeights();
};

#endif
