#ifndef LEGENDRE_H
#define LEGENDRE_H

#include "polynomial.h"
#include "mathtools.h"

#include <boost/math/special_functions/legendre.hpp>

class Legendre : public Polynomial{
private:
    void ComputeNodes(int degree);
public:
    Legendre(int degree);

    virtual double Evaluate(int m, double x);
    virtual blaze::DynamicVector<double> GetNodes();
    virtual blaze::DynamicVector<double> GetWeights();
};

#endif
