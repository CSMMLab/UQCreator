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
    virtual MatVec InitLambda( const MatVec& u );
    virtual double ComputeDt( Vector& u, double dx ) const;
};

#endif    // EULER2D_H
