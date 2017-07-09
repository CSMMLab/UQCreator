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
    double _dx,_dt,_a,_b,_uL,_uR;
    vector numFlux(vector lambda1,vector lambda2);
    std::vector<vector> SetupIC();
    double IC(double x,double xi, double uL, double uR);
public:
    MomentSolver(Problem* problem);
    void Solve();
};

#endif // MOMENTSOLVER_H
