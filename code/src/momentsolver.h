#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <blaze/math/DynamicVector.h>
#include <gnuplot-iostream.h>
#include <chrono>

#include "closure.h"
#include "problem.h"
#include "legendre.h"

class MomentSolver
{
    Legendre* _quad;
    Closure* _closure;
    Mesh* _mesh;
    blaze::DynamicVector<double> _x;
    std::vector<blaze::DynamicVector<double> > _lambda;
    Problem* _problem;
    double _dx,_dt,_a,_b,_uL,_uR,_tEnd;
    int _nTimeSteps,_nCells,_nMoments;
    blaze::DynamicVector<double> numFlux(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2);
    std::vector<blaze::DynamicVector<double>> SetupIC();
    double IC(double x,double xi);
    blaze::DynamicVector<double> CalculateMoments(const blaze::DynamicVector<double>& lambda);
    blaze::DynamicVector<double> EvalLambda(const blaze::DynamicVector<double>& lambda, const blaze::DynamicVector<double>& xi);
public:
    MomentSolver(Problem* problem);
    void Solve();
    void Plot();
    void Print();
    void PlotFixedXi();
};

#endif // MOMENTSOLVER_H

