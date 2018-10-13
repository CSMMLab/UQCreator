#ifndef EXPLICITEULER_H
#define EXPLICITEULER_H

#include "timesolver.h"

class ExplicitEuler : public TimeSolver
{
  private:
  public:
    ExplicitEuler() = delete;
    ExplicitEuler( Problem* problem );
    virtual void Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                          std::vector<Matrix>& uNew,
                          std::vector<Matrix>& u,
                          std::vector<Matrix>& lambda );
};

#endif    // EXPLICITEULER_H
