#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <iostream>

class TimeSolver
{
private:

protected:

public:
    TimeSolver();
    virtual ~TimeSolver();
    static TimeSolver* Create(std::string inputFile);
    virtual double Solve(const blaze::DynamicVector<double>& u, const blaze::DynamicVector<double>& flux1, const blaze::DynamicVector<double>& flux2) = 0;
};

#endif // TIMEDISCRETIZATION_H
