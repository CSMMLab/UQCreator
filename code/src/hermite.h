#ifndef HERMITE_H
#define HERMITE_H

#include "polynomial.h"
#include "mathtools.h"

#include <boost/math/special_functions/hermite.hpp>

class Hermite : public Polynomial{
private:
    virtual void Compute();
public:
    Hermite(int degree);

    virtual double Evaluate(int m,double x);
    virtual const blaze::DynamicVector<double>& GetNodes();
    virtual const blaze::DynamicVector<double>& GetWeights();
};

#endif
