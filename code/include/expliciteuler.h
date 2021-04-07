#ifndef EXPLICITEULER_H
#define EXPLICITEULER_H

#include "problems/euler2d.h"
#include "timesolver.h"

class ExplicitEuler : public TimeSolver
{
  private:
    ExplicitEuler() = delete;
    std::vector<Cell*> _cells;

  public:
    ExplicitEuler( Settings* settings, Mesh* mesh, Problem* problem );
    virtual void Advance( std::function<void( Matrix&, const Matrix&, unsigned )> const& fluxFunc,
                          MatTens& uNew,
                          MatTens& u,
                          MatTens& uQ,
                          double dt,
                          const VectorU& refLevel );
};

#endif    // EXPLICITEULER_H
