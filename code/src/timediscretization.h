#ifndef TIMEDISCRETIZATION_H
#define TIMEDISCRETIZATION_H

#include <blaze/Blaze.h>
#include <list>

typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;

class TimeDiscretization
{
    vector _alpha;
    vector _beta;
    std::list<vector> _uSteps;
public:
    TimeDiscretization();
    vector Solve(vector u, vector Lambda);
};

#endif // TIMEDISCRETIZATION_H
