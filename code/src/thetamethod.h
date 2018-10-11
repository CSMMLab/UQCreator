#ifndef THETAMETHOD_H
#define THETAMETHOD_H

#include "timesolver.h"

class ThetaMethod : public TimeSolver
{
  private:
    double _theta;

  public:
    ThetaMethod() = delete;
    ThetaMethod( Problem* problem, double theta );
    void Advance( std::function<Matrix( const Matrix&, const Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                  std::vector<Matrix>& uNew,
                  std::vector<Matrix>& u,
                  std::vector<Matrix>& lambda );
};

#endif    // THETAMETHOD_H
