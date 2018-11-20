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
    Legendre( unsigned degree );
    virtual ~Legendre() {}

    virtual double Evaluate( unsigned m, double x );
    virtual const Vector& GetNodes();
    virtual const Vector& GetWeights();
    virtual double fXi( const double xi ) const;
};

#endif
