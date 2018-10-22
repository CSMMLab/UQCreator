#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <cpptoml.h>
#include <functional>
#include <iostream>

#include "mesh.h"
#include "settings.h"

class TimeSolver
{
  private:
  protected:
    const Settings* _settings;
    const Mesh* _mesh;
    double _CFL;
    double _dt;
    double _dx;
    double _tEnd;
    unsigned _nTimeSteps;

    TimeSolver() = delete;

  public:
    TimeSolver( Settings* settings, Mesh* mesh );
    virtual ~TimeSolver();
    static TimeSolver* Create( Settings* settings, Mesh* mesh );
    virtual void Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                          std::vector<Matrix>& uNew,
                          std::vector<Matrix>& u,
                          std::vector<Matrix>& uQ ) = 0;
    double GetTimeStepSize();
    double GetNTimeSteps();
};

#endif    // TIMEDISCRETIZATION_H
