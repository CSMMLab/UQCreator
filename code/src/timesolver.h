#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <functional>
#include <iostream>

#include "closure.h"
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
    TimeSolver() = delete;
    TimeSolver( Problem* problem );
    virtual ~TimeSolver();
    static TimeSolver* Create( Problem* problem, Closure* closure );
    virtual void Advance( std::function<Matrix( const Matrix&,
                                                                      const Matrix&,
                                                                      const Matrix&,
                                                                      const Matrix& )> const& fluxFunc,
                          std::vector<Matrix>& uNew,
                          std::vector<Matrix>& u,
                          std::vector<Matrix>& lambda ) = 0;
    double GetTimeStepSize();
    double GetNTimeSteps();
};

#endif    // TIMEDISCRETIZATION_H
