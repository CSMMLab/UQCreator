#ifndef BURGERS_H
#define BURGERS_H

#include <fstream>
#include <iostream>

#include "gnuplot-iostream.h"
#include "problem.h"

class Burgers : public Problem
{
  private:
    Vector _u;
    Vector _x;
    double _dx;
    double _dt;
    unsigned _nCells;
    unsigned _nTimeSteps;
    double H( double u, double v, double w );
    double F( double u );
    Matrix F( const Matrix& u );
    double IC( double x, double uL, double uR );

    Burgers() {}

  public:
    Burgers( std::string inputFile );
    virtual void Solve();
    virtual void Plot( Vector& x, Vector& u );
    virtual void Print();
    virtual void WriteToFile( std::string filename, int filetype ) const;
    double G( double u, double v );
    Matrix G( const Matrix& u, const Matrix& v );
    virtual double ExactSolution( double t, double x, double xi );
};

#endif    // BURGERS_H
