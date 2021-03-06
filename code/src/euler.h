#ifndef EULER_H
#define EULER_H

#include "problem.h"

class Euler : public Problem
{
  private:
    double _gamma;
    double SodFunction( double P, double rho_l, double P_l, double rho_r, double P_r ) const;
    double Bisection( double PA, double PB, double rho_l, double P_l, double rho_r, double P_r ) const;

  public:
    Euler( Settings* settings );
    virtual ~Euler();
    virtual void Solve();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
};

#endif    // EULER_H
