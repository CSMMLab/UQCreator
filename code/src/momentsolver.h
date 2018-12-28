#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <chrono>
#include <omp.h>
#include <spdlog/spdlog.h>

#include "closure.h"
#include "mesh.h"
#include "polynomial.h"
#include "problem.h"
#include "settings.h"
#include "timesolver.h"
#include "typedefs.h"

class MomentSolver
{
  private:
    Settings* _settings;
    Closure* _closure;
    Mesh* _mesh;
    TimeSolver* _time;
    MatVec _lambda;
    Problem* _problem;
    double _dt, _tEnd;
    unsigned _nCells;         // number of spatial cells
    unsigned _nMoments;       // number of moments in one uncertain dimension
    unsigned _nStates;        // number of states of the original system
    unsigned _nQuadPoints;    // number of moments in one uncertain dimension
    unsigned _nQTotal;        // total number of quad points
    unsigned _nTotal;         // total number of moments
    std::shared_ptr<spdlog::logger> _log;

    void numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n );
    void SetupIC( MatVec& out );
    MatVec SetupIC();
    void CalculateMoments( MatVec& out, const MatVec& lambda );
    Vector EvalLambda( const Vector& lambda, const Vector& xi );
    void Plot( double time, unsigned nSteps );

  public:
    MomentSolver( Settings* settings, Mesh* mesh, Problem* problem );
    ~MomentSolver();
    void Solve();
};

#endif    // MOMENTSOLVER_H
