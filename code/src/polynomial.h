#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

class Polynomial
{
protected:
    int _degree;
    blaze::DynamicVector<double> _nodes;
    blaze::DynamicVector<double> _weights;

    void Sort();
    virtual void Compute()=0;
public:
    Polynomial(int degree);
    virtual double Evaluate(int m, double x)=0;
    virtual const blaze::DynamicVector<double>& GetNodes()=0;
    virtual const blaze::DynamicVector<double>& GetWeights()=0;
};

#endif // POLYNOMIAL_H
