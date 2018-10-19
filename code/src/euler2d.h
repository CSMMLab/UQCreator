#ifndef EULER2D_H
#define EULER2D_H

#include "problem.h"

class Euler2D : public Problem
{
  private:
    double _gamma;

  public:
    Euler2D( const Settings* settings );
    virtual ~Euler2D();
    virtual void Solve();
    virtual Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    virtual double ExactSolution( double t, double x, double xi );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    double GetGamma() const;
};

#endif    // EULER2D_H
