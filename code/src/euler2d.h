#ifndef EULER2D_H
#define EULER2D_H

#include "problem.h"

class Euler2D : public Problem
{
  private:
    double _gamma;

  public:
    Euler2D( Settings* settings );
    virtual ~Euler2D();
    virtual void Solve();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif    // EULER2D_H
