#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include <blaze/math/DynamicVector.h>

#include "closure.h"
#include "problem.h"
#include "hermite.h"
#include "gnuplot-iostream.h"

class MomentSolver
{
    Polynomial* _quad;
    Closure* _closure;
    Mesh* _mesh;
    blaze::DynamicVector<double> _x;
    std::vector<blaze::DynamicVector<double> > _lambda;
    Problem* _problem;
    double _dx,_dt,_a,_b,_uL,_uR,_tEnd;
    int _nTimeSteps,_nCells,_nMoments;
    blaze::DynamicVector<double> numFlux(blaze::DynamicVector<double> lambda1,blaze::DynamicVector<double> lambda2);
    std::vector<blaze::DynamicVector<double>> SetupIC();
    double IC(double x,double xi);
    blaze::DynamicVector<double> CalculateMoments(blaze::DynamicVector<double> lambda);
    blaze::DynamicVector<double> EvalLambda(blaze::DynamicVector<double> lambda, blaze::DynamicVector<double> xi);
public:
    MomentSolver(Problem* problem);
    void Solve();
    void Plot();
    void Print();
};

#endif // MOMENTSOLVER_H

