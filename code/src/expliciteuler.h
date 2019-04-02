#ifndef EXPLICITEULER_H
#define EXPLICITEULER_H

#include "euler2d.h"
#include "timesolver.h"

class ExplicitEuler : public TimeSolver
{
  private:
    ExplicitEuler() = delete;
    std::vector<Cell*> _cells;
    Matrix _ghostCell;

  public:
    ExplicitEuler( Settings* settings, Mesh* mesh );
    virtual void Advance( std::function<void( Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector&, unsigned )> const& fluxFunc,
                          MatVec& uNew,
                          MatVec& u,
                          MatVec& uQ,
                          double dt,
                          const std::vector<unsigned>& refLevel );
};

#endif    // EXPLICITEULER_H
