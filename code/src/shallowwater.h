#ifndef SHALLOWWATER_H
#define SHALLOWWATER_H

#include "problem.h"

class ShallowWater : public Problem
{
  private:
    double _g;

  public:
    ShallowWater( Settings* settings );
    virtual ~ShallowWater();
    virtual void Solve();
    inline Vector G( const Vector& q_l, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif    // SHALLOWWATER_H
