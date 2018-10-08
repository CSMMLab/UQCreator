#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <blaze/math/DynamicVector.h>
#include <chrono>

#include "closure.h"
#include "legendre.h"
#include "limiter.h"
#include "minmod.h"
#include "nolimiter.h"
#include "plotengine.h"
#include "problem.h"
#include "thetamethod.h"
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
    double _dx, _dt, _a, _b, _uL, _uR, _tEnd;
    unsigned _nTimeSteps, _nCells, _nMoments, _nStates;
    Matrix numFlux( const Matrix& lambda0, const Matrix& lambda1, const Matrix& lambda2, const Matrix& lambda3 );
    std::vector<Matrix> SetupIC();
    double IC( double x, double xi );
    Matrix CalculateMoments( const Matrix& lambda );
    Vector EvalLambda( const Vector& lambda, const Vector& xi );

  public:
    MomentSolver( Problem* problem );
    ~MomentSolver();
    void Solve();
    void Print();
    Result1D GetPlotData1D();
    Result1D GetPlotData1DFixedXi();
    Result1D GetPlotData1DExpectedValue();
};

#endif    // MOMENTSOLVER_H
