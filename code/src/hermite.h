#ifndef HERMITE_H
#define HERMITE_H

#include "polynomial.h"
#include "mathtools.h"

class Hermite : public Polynomial{
private:
    void computeNodes(int degree);
public:
    Hermite(int degree);

    virtual double evaluate(int m,double x);
    virtual vector getNodes();
    virtual vector getWeights();
};

#endif
