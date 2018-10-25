#ifndef EXPLICITEULER_H
#define EXPLICITEULER_H

#include "euler2d.h"
#include "timesolver.h"

class ExplicitEuler : public TimeSolver
{
  private:
    ExplicitEuler() = delete;

  public:
    ExplicitEuler( Settings* settings, Mesh* mesh );
    virtual void Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                          std::vector<Matrix>& uNew,
                          std::vector<Matrix>& u,
                          std::vector<Matrix>& uQ );
};

#endif    // EXPLICITEULER_H
