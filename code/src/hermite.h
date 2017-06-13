#ifndef HERMITE_H
#define HERMITE_H

#include "mathtools.h"

class Hermite : public Quadrature{
private:
    void computeNodes(int degree);
public:
    virtual double evaluate();
    virtual vector getRoots();
    virtual vector getWeights();
};

#endif
