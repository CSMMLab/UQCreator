#ifndef KINETICEQUATION_H
#define KINETICEQUATION_H

#include "problem.h"
class QuadratureGrid;

class KineticEquation : public Problem
{
  private:
    std::vector<Vector> _xiQuad;
    QuadratureGrid* _grid;

  public:
    KineticEquation( Settings* settings );
    virtual ~KineticEquation();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u, double mu );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
};

#endif    // KINETICEQUATION_H
