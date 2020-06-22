#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <cpptoml.h>
#include <functional>
#include <iostream>
#include <omp.h>
#include <spdlog/spdlog.h>

#include "mesh.h"
#include "problem.h"
#include "settings.h"

class TimeSolver
{
  private:
  protected:
    const Settings* _settings;
    const Mesh* _mesh;
    Problem* _problem;
    double _CFL;
    double _dt;
    double _dx;
    double _tEnd;
    unsigned _nTimeSteps;
    MatVec _flux;    // flux over cell edges

    TimeSolver() = delete;

  public:
    TimeSolver( Settings* settings, Mesh* mesh, Problem* problem );
    virtual ~TimeSolver();
    static TimeSolver* Create( Settings* settings, Mesh* mesh, Problem* problem );
    virtual void Advance( std::function<void( Matrix&, const Matrix&, unsigned )> const& fluxFunc,
                          MatTens& uNew,
                          MatTens& u,
                          MatTens& uQ,
                          double dt,
                          const VectorU& refLevel ) = 0;
    double GetTimeStepSize();
    double GetNTimeSteps();
};

#endif    // TIMEDISCRETIZATION_H
