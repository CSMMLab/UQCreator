#ifndef EULER2D_H
#define EULER2D_H

#include "problem.h"

enum ICType { I_NOZZLE, I_NOZZLE_SOD, I_NACA, I_NACA_HIGHMACH };

class Euler2D : public Problem
{
  private:
    double _gamma;
    ICType _problemType;

  public:
    Euler2D( Settings* settings );
    virtual ~Euler2D();
    virtual void Solve();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u ) const;
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Tensor& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
    virtual Matrix BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const;
};

#endif    // EULER2D_H