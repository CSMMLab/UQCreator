#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "quadrature.h"
#include "mathtools.h"

class Legendre : private Quadrature{
private:
    void computeNodes(int degree);
public:
    Legendre(int degree);

    virtual double evaluate();
    virtual vector getNodes();
    virtual vector getWeights();
};

#endif
