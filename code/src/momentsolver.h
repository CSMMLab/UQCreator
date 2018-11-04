#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <chrono>
#include <omp.h>

#include "closure.h"
#include "legendre.h"
#include "mesh.h"
#include "problem.h"
#include "settings.h"
#include "timesolver.h"
#include "typedefs.h"

class MomentSolver
{
  private:
    Settings* _settings;
    Legendre* _quad;
    Closure* _closure;
    Mesh* _mesh;
    TimeSolver* _time;
    MatVec _lambda;
    Problem* _problem;
    double _dt, _tEnd;
    unsigned _nCells, _nMoments, _nStates, _nQuadPoints;

    void numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n );
    void SetupIC( MatVec& out );
    Vector IC( Vector x, double xi );
    void CalculateMoments( MatVec& out, const MatVec& lambda );
    Vector EvalLambda( const Vector& lambda, const Vector& xi );
    void Plot( double time, unsigned nSteps );

  public:
    MomentSolver( Settings* settings, Mesh* mesh, Problem* problem );
    ~MomentSolver();
    void Solve();
};

#endif    // MOMENTSOLVER_H
