#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/Blaze.h>

#include "problem.h"
#include "polynomial.h"

typedef blaze::DynamicVector<double> vector;

class Quadrature{
protected:
    double _value;
    vector _nodes;
    vector _weights;

    Problem* _problem;
    Polynomial* _polynomial;


public:
    Quadrature(Problem* p);
    double evaluate();
    vector getNodes();
    vector getWeights();
private:
    Quadrature(){}
};

#endif // QUADRATURE_H
