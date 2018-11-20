#ifndef HERMITE_H
#define HERMITE_H

#include "mathtools.h"
#include "polynomial.h"

#include <boost/math/special_functions/hermite.hpp>

class Hermite : public Polynomial
{
  private:
    virtual void Compute();

  public:
    Hermite( unsigned degree );

    virtual double Evaluate( unsigned m, double x );
    virtual const Vector& GetNodes();
    virtual const Vector& GetWeights();
    virtual double fXi( const double xi ) const;
    virtual double L2NormSquare( unsigned i ) const;
    unsigned factorial( unsigned n ) const;
};

#endif
