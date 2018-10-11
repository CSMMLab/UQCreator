#ifndef SSPMULTISTEP_H
#define SSPMULTISTEP_H

#include "heun.h"
#include "timesolver.h"

class SSPMultiStep : public TimeSolver
{
    std::vector<Matrix> _u1Step, _u2Step, _u3Step;
    Heun* _heun;
    unsigned int _counter;

  public:
    SSPMultiStep() = delete;
    SSPMultiStep( Problem* problem, Closure* closure );
    virtual ~SSPMultiStep();
    void Advance( std::function<Matrix( const Matrix&, const Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                  std::vector<Matrix>& uNew,
                  std::vector<Matrix>& u,
                  std::vector<Matrix>& lambda ) {}
};

#endif    // SSPMULTISTEP_H
