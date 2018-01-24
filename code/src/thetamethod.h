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
    virtual void Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                          std::vector<blaze::DynamicMatrix<double>>& uNew,
                          const std::vector<blaze::DynamicMatrix<double>>& u,
                          const std::vector<blaze::DynamicMatrix<double>>& lambda );
};

#endif    // THETAMETHOD_H
