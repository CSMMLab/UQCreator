#ifndef SSPMULTISTEP_H
#define SSPMULTISTEP_H

#include "timesolver.h"
#include "heun.h"

class SSPMultiStep : public TimeSolver
{
    std::vector<blaze::DynamicMatrix<double>> _u1Step, _u2Step, _u3Step;
    Heun* _heun;
    unsigned int _counter;
public:
    SSPMultiStep() = delete;
    SSPMultiStep( Problem * problem, Closure* closure);
    virtual ~SSPMultiStep();
    virtual void Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>&,
                                                                      const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                          std::vector<blaze::DynamicMatrix<double>>& uNew,
                          std::vector<blaze::DynamicMatrix<double>>& u,
                          std::vector<blaze::DynamicMatrix<double>>& lambda );
};

#endif // SSPMULTISTEP_H
