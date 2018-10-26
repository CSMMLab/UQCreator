#ifndef BURGERS_H
#define BURGERS_H

#include <fstream>
#include <iostream>

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
    Vector F( double u );
    Matrix F( const Matrix& u );
    double IC( double x, double uL, double uR );

    Burgers() {}

  public:
    Burgers( Settings* settings );
    Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    virtual double ExactSolution( double t, double x, double xi );
    virtual double ComputeDt( Vector& u, double dx ) const;
};

#endif    // BURGERS_H
