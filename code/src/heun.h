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
    virtual void Advance( std::function<Matrix( const Matrix&,
                                                                      const Matrix&,
                                                                      const Matrix&,
                                                                      const Matrix& )> const& fluxFunc,
                          std::vector<Matrix>& uNew,
                          std::vector<Matrix>& u,
                          std::vector<Matrix>& lambda );
};

#endif // HEUN_H
