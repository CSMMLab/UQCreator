#ifndef SHALLOWWATER2D_H
#define SHALLOWWATER2D_H

#include "problem.h"

class ShallowWater2D : public Problem
{
  private:
    double _g;    // gravitational constant

  public:
    ShallowWater2D( Settings* settings );
    virtual ~ShallowWater2D();
    virtual void Solve();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif    // SHALLOWWATER2D_H
