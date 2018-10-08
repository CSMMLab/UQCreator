#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

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
    virtual double Evaluate( unsigned m, double x ) = 0;
    virtual const Vector& GetNodes()                = 0;
    virtual const Vector& GetWeights()              = 0;
};

#endif    // POLYNOMIAL_H
