#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/math/DynamicVector.h>

#include "polynomial.h"
#include "problem.h"

class Quadrature
{
  protected:
    double _value;
    Vector _nodes;
    Vector _weights;

    Problem* _problem;
    Polynomial* _polynomial;

  public:
    Quadrature( Problem* p );
    double Evaluate();
    Vector GetNodes();
    Vector GetWeights();

  private:
    Quadrature() {}
};

#endif    // QUADRATURE_H
