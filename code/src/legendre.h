#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "polynomial.h"
#include "mathtools.h"

class Legendre : public Polynomial{
private:
    void ComputeNodes(int degree);
public:
    Legendre(int degree);

    virtual double Evaluate(int m, double x);
    virtual vector GetNodes();
    virtual vector GetWeights();
};

#endif
