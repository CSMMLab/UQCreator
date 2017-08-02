#ifndef BURGERS_H
#define BURGERS_H

#include <blaze/Blaze.h>
#include <iostream>
#include <fstream>

#include "problem.h"
#include "gnuplot-iostream.h"

typedef blaze::DynamicVector<double> vector;

class Burgers : public Problem {
private:
    vector _u;
    vector _x;
    double _dx;
    double _dt;
    int _nCells;
    int _nTimeSteps;
    double H(double u, double v, double w);
    double F(double u);
    double IC(double x, double uL, double uR);

    //Burgers(){}

public:
    Burgers(std::string inputFile);
    Burgers(){}
    virtual void Solve();
    virtual void Plot() const;
    virtual void Print() const;
    virtual void WriteToFile(std::string filename, int filetype) const;
    double G(double u, double v);
};

#endif // BURGERS_H