#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <algorithm>
#include <numeric>

#include "settings.h"
#include "typedefs.h"

class Polynomial
{
  protected:
    unsigned _degree;
    Vector _nodes;
    Vector _weights;

    void Sort();
    virtual void Compute() = 0;

  public:
    Polynomial( unsigned degree );
    virtual ~Polynomial() {}
    static Polynomial* Create( Settings* settings, unsigned order );
    virtual double Evaluate( unsigned m, double x ) = 0;
    virtual const Vector& GetNodes()                = 0;
    virtual const Vector& GetWeights()              = 0;
    virtual double fXi( const double xi ) const     = 0;
};

#endif    // POLYNOMIAL_H
