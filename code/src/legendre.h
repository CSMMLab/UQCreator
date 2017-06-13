#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "mathtools.h"

class Legendre : public Quadrature{
private:
    void computeNodes(int degree);
public:
    virtual double evaluate();
    virtual vector getRoots();
    virtual vector getWeights();
};

#endif
