#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/math/DynamicVector.h>

#include "problem.h"
#include "polynomial.h"

class Quadrature{
protected:
    double _value;
    blaze::DynamicVector<double> _nodes;
    blaze::DynamicVector<double> _weights;

    Problem* _problem;
    Polynomial* _polynomial;


public:
    Quadrature(Problem* p);
    double Evaluate();
    blaze::DynamicVector<double> GetNodes();
    blaze::DynamicVector<double> GetWeights();
private:
    Quadrature(){}
};

#endif // QUADRATURE_H
