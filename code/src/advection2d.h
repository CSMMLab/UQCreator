#ifndef ADVECTION2D_H
#define ADVECTION2D_H

#include <fstream>
#include <iostream>

#include "problem.h"

class Advection2D : public Problem
{
  private:
    // parameters for initial condition
    double _uL;
    double _uR;
    double _x0;
    double _x1;

    Vector _omega; // vector with advection direction

    Vector F( double u );
    Matrix F( const Matrix& u );

  public:
    Advection2D( Settings* settings );
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    virtual double ComputeDt(const Tensor &u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
};

#endif // ADVECTION2D_H
