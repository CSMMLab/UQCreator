#ifndef TIMEDISCRETIZATION_H
#define TIMEDISCRETIZATION_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <list>

class TimeDiscretization
{
    blaze::DynamicVector<double> _alpha;
    blaze::DynamicVector<double> _beta;
    std::list<blaze::DynamicVector<double>> _uSteps;
public:
    TimeDiscretization();
    blaze::DynamicVector<double> Solve(blaze::DynamicVector<double> u, blaze::DynamicVector<double> Lambda);
};

#endif // TIMEDISCRETIZATION_H
