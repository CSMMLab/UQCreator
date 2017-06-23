#ifndef BURGERSSOLVER_H
#define BURGERSSOLVER_H

#include <iostream>
#include <fstream>

#include "gnuplot-iostream.h"


class BurgersSolver
{
    double* _u;
    double* _x;
    int _nTimeSteps;
    int _nCells;
    double _dt, _dx;
    double H(double u, double v, double w);
    double g(double u, double v);
    double f(double u);
    double IC(double x, double uL, double uR);
    double _a,_b;
public:
    BurgersSolver(int nCells, double tEnd, double cfl, double a, double b, double uL, double uR);
    void Solve();
    void Print() const;
    void Plot() const;
};

#endif // BURGERSSOLVER_H
