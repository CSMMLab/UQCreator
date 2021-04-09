#ifndef EXTERNALPROBLEM_H
#define EXTERNALPROBLEM

#include "problem.h"

class ExternalProblem : public Problem
{
  private:

  public:
    ExternalProblem( Settings* settings );
    virtual ~ExternalProblem();
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

#endif // EXTERNALPROBLEM
