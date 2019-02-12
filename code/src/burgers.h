#ifndef BURGERS_H
#define BURGERS_H

#include <fstream>
#include <iostream>

#include "problem.h"

class Burgers : public Problem
{
  private:
    // parameters for initial condition
    double _uL;
    double _uR;
    double _x0;
    double _x1;

    Vector F( double u );
    Matrix F( const Matrix& u );

    Burgers() {}

  public:
    Burgers( Settings* settings );
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    virtual double ComputeDt( Vector& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
    Vector ExactSolution( double t, const Vector& x, const Vector& xi ) const;
};

#endif    // BURGERS_H
