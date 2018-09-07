#ifndef HEUN_H
#define HEUN_H

#include "timesolver.h"
#include "closure.h"


class Heun : public TimeSolver
{
    Closure* _closure;
public:
    Heun() = delete;
    Heun( Problem* problem, Closure* closure );
    virtual ~Heun();
    virtual void Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                          std::vector<blaze::DynamicMatrix<double>>& uNew,
                          std::vector<blaze::DynamicMatrix<double>>& u,
                          std::vector<blaze::DynamicMatrix<double>>& lambda );
};

#endif // HEUN_H
