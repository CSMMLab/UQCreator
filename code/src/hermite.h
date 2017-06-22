#ifndef HERMITE_H
#define HERMITE_H

#include "quadrature.h"
#include "mathtools.h"

class Hermite : public Quadrature{
private:
    void computeNodes(int degree);
public:
    Hermite(int degree);

    virtual double evaluate();
    virtual vector getNodes();
    virtual vector getWeights();
};

#endif
