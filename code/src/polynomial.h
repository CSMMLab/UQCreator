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
    virtual double Evaluate(int m, double x)=0;
    virtual vector GetNodes()=0;
    virtual vector GetWeights()=0;
};

#endif // POLYNOMIAL_H
