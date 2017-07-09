#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <blaze/Blaze.h>

typedef blaze::DynamicVector<double> vector;

class Polynomial
{
protected:
    int _degree;
    vector _nodes;
    vector _weights;
public:
    Polynomial(int degree);
    virtual double evaluate(int m, double x)=0;
    virtual vector getNodes()=0;
    virtual vector getWeights()=0;
};

#endif // POLYNOMIAL_H
