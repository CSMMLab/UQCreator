#ifndef EULER_H
#define EULER_H

#include "problem.h"

class Euler : public Problem
{
  private:
    double _gamma;

  public:
    Euler( Settings* settings );
    virtual ~Euler();
    virtual void Solve();
    virtual Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    virtual double ExactSolution( double t, double x, double xi );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( Vector& u, double dx ) const;
};

#endif    // EULER_H
