#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H

#include "problem.h"

class NavierStokes : public Problem
{
private:
  double _gamma;
public:
  NavierStokes( Settings* settings );
  virtual ~NavierStokes();
  inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
  virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
  //virtual Matrix Source( const Matrix& uQ ) const;
  Matrix F( const Vector& u );
  //Matrix F( const Matrix& u );
  //virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
  virtual Vector IC( const Vector& x, const Vector& xi );
  virtual Vector LoadIC( const Vector& x, const Vector& xi ) {};
  //virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
};

#endif // NAVIERSTOKES_H
