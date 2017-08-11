#ifndef BURGERS_H
#define BURGERS_H

#include <iostream>
#include <fstream>

#include "problem.h"
#include "gnuplot-iostream.h"

class Burgers : public Problem {
private:
    blaze::DynamicVector<double> _u;
    blaze::DynamicVector<double> _x;
    double _dx;
    double _dt;
    int _nCells;
    int _nTimeSteps;
    double H(double u, double v, double w);
    double F(double u);
    double IC(double x, double uL, double uR);

    Burgers(){}
public:
    Burgers(std::string inputFile);
    virtual void Solve();
    virtual void Plot(blaze::DynamicVector<double>& x, blaze::DynamicVector<double>& u) const;
    virtual void Print() const;
    virtual void WriteToFile(std::string filename, int filetype) const;
    double G(double u, double v);
};

#endif // BURGERS_H
