#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <blaze/Blaze.h>

#include "problem.h"

typedef blaze::DynamicVector<double> vector;

class Quadrature{
protected:
    double _value;
    vector _nodes;
    vector _weights;

    Problem* _problem;

public:
    Quadrature(){} // move to private
    Quadrature(Problem* p);
    virtual double evaluate() = 0;
    virtual vector getNodes() = 0;
    virtual vector getWeights() = 0;
};

#endif // QUADRATURE_H
