#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <blaze/math/DynamicVector.h>
#include <chrono>

#include "closure.h"
#include "expliciteuler.h"
#include "legendre.h"
#include "limiter.h"
#include "minmod.h"
#include "nolimiter.h"
#include "omp.h"
#include "plotengine.h"
#include "problem.h"
#include "typedefs.h"

class MomentSolver
{
  private:
    Legendre* _quad;
    Closure* _closure;
    Mesh* _mesh;
    TimeSolver* _time;
    Limiter* _limiter;
    Vector _x;
    std::vector<Matrix> _lambda;
    Problem* _problem;
    PlotEngine* _plotEngine;
    double _dx, _dt, _a, _b, _uL, _uR, _tEnd;
    unsigned _nTimeSteps, _nCells, _nMoments, _nStates;
    Matrix numFlux( const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n );
    std::vector<Matrix> SetupIC();
    Vector IC( Vector x, double xi );
    Matrix CalculateMoments( const Matrix& lambda );
    Vector EvalLambda( const Vector& lambda, const Vector& xi );
    void Plot( double time );

  public:
    MomentSolver( Problem* problem );
    ~MomentSolver();
    void Solve();
};

#endif    // MOMENTSOLVER_H
