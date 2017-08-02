#ifndef MOMENTSOLVER_H
#define MOMENTSOLVER_H

#include "burgers.h"
#include "closure.h"
#include "problem.h"
#include "hermite.h"

typedef blaze::DynamicVector<double> vector;

class MomentSolver
{
    Polynomial* _quad;
    Closure* _closure;
    vector _x;
    Problem* _problem;    
    Burgers* _origSolver;
    double _dx,_dt,_a,_b,_uL,_uR,_nCells,_nMoments,_tEnd;
    vector numFlux(vector lambda1,vector lambda2);
    std::vector<vector> SetupIC();
    double IC(double x,double xi);
    vector CalculateMoments(vector lambda);
public:
    MomentSolver(Problem* problem);
    void Solve();
};

#endif // MOMENTSOLVER_H
