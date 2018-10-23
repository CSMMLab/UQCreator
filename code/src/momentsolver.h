#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <chrono>
#include <omp.h>

#include "closure.h"
#include "expliciteuler.h"
#include "legendre.h"
#include "limiter.h"
#include "mesh.h"
#include "minmod.h"
#include "nolimiter.h"
#include "plotengine.h"
#include "problem.h"
#include "settings.h"
#include "typedefs.h"

class MomentSolver
{
  private:
    Settings* _settings;
    Legendre* _quad;
    Closure* _closure;
    Mesh* _mesh;
    TimeSolver* _time;
    Limiter* _limiter;
    std::vector<Matrix> _lambda;
    Problem* _problem;
    PlotEngine* _plotEngine;
    double _dt, _tEnd;
    unsigned _nCells, _nMoments, _nStates, _nQuadPoints;
    Matrix numFlux( const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n );
    std::vector<Matrix> SetupIC();
    Vector IC( Vector x, double xi );
    Matrix CalculateMoments( const Matrix& lambda );
    Vector EvalLambda( const Vector& lambda, const Vector& xi );
    void Plot( double time, unsigned nSteps );

  public:
    MomentSolver( Settings* settings, Mesh* mesh, Problem* problem );
    ~MomentSolver();
    void Solve();
};

#endif    // MOMENTSOLVER_H
