#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <functional>
#include <iostream>

#include "problem.h"

class TimeSolver
{
  private:
  protected:
    Problem* _problem;
    double _CFL;
    double _dt;
    double _dx;
    double _tEnd;
    int _nTimeSteps;

  public:
    TimeSolver( Problem* problem );
    virtual ~TimeSolver();
    static TimeSolver* Create( Problem* problem );
    virtual void Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                          std::vector<blaze::DynamicMatrix<double>>& uNew,
                          const std::vector<blaze::DynamicMatrix<double>>& u,
                          const std::vector<blaze::DynamicMatrix<double>>& lambda ) = 0;
    double GetTimeStepSize();
    double GetNTimeSteps();
};

#endif    // TIMEDISCRETIZATION_H
