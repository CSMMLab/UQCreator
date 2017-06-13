#ifndef QUADRATURE_H
#define QUADRATURE_H

typedef blaze::DynamicVector<double> vector;

class Quadrature{
private:
    double _value;
    vector _roots;
    vector _weights;

    Problem* _problem;

    Quadrature(){}
public:
    Quadrature(Problem* p);
    virtual double evaluate() = 0;
    virtual vector getRoots() = 0;
    virtual vector getWeights() = 0;
};

#endif // QUADRATURE_H
