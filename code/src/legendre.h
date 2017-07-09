#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "polynomial.h"
#include "mathtools.h"

class Legendre : public Polynomial{
private:
    void computeNodes(int degree);
public:
    Legendre(int degree);

    virtual double evaluate(int m, double x);
    virtual vector getNodes();
    virtual vector getWeights();
};

#endif
