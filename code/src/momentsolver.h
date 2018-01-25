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

class MomentSolver
{
  private:
    Legendre* _quad;
    Closure* _closure;
    Mesh* _mesh;
    TimeSolver* _time;
    Limiter* _limiter;
    blaze::DynamicVector<double> _x;
    std::vector<blaze::DynamicMatrix<double>> _lambda;
    Problem* _problem;
    PlotEngine* _plot;
    double _dx, _dt, _a, _b, _uL, _uR, _tEnd;
    int _nTimeSteps, _nCells, _nMoments, _nStates;
    blaze::DynamicMatrix<double> numFlux( const blaze::DynamicMatrix<double>& lambda0,
                                          const blaze::DynamicMatrix<double>& lambda1,
                                          const blaze::DynamicMatrix<double>& lambda2,
                                          const blaze::DynamicMatrix<double>& lambda3 );
    std::vector<blaze::DynamicMatrix<double>> SetupIC();
    double IC( double x, double xi );
    blaze::DynamicMatrix<double> CalculateMoments( const blaze::DynamicMatrix<double>& lambda );
    blaze::DynamicVector<double> EvalLambda( const blaze::DynamicVector<double>& lambda, const blaze::DynamicVector<double>& xi );

  public:
    MomentSolver( Problem* problem );
    ~MomentSolver();
    void Solve();
    void Plot();
    void Print();
    void PlotFixedXi();
};

#endif    // MOMENTSOLVER_H
