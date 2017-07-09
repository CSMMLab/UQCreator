#ifndef HERMITE_H
#define HERMITE_H

#include "polynomial.h"
#include "mathtools.h"

class Hermite : public Polynomial{
private:
    void ComputeNodes(int degree);
public:
    Hermite(int degree);

    virtual double Evaluate(int m,double x);
    virtual vector GetNodes();
    virtual vector GetWeights();
};

#endif
