#ifndef THETAMETHOD_H
#define THETAMETHOD_H

#include "timesolver.h"

class ThetaMethod : public TimeSolver
{
private:
    double _theta;
public:
    ThetaMethod() = delete;
    ThetaMethod(double theta);
    double Solve(const blaze::DynamicVector<double>& u, const blaze::DynamicVector<double>& flux1, const blaze::DynamicVector<double>& flux2);
};

#endif // THETAMETHOD_H
