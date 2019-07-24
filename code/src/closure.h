#ifndef CLOSURE_H
#define CLOSURE_H

#include <spdlog/spdlog.h>
#include <vector>

#include "legendre.h"
#include "settings.h"
#include "typedefs.h"

#include "quadraturegrid.h"

class Closure
{
  private:
    Closure() = delete;

  protected:
    Settings* _settings;
    std::vector<Polynomial*> _quad;
    Matrix _phiTildeWf;    // stores scaled basis functions evaluated at quadrature points times weight and pdf
    double _alpha;         // step size for Newton
    unsigned _nMoments;
    unsigned _nQuadPoints;
    unsigned _nStates;
    unsigned _numDimXi;
    unsigned _nQTotal;
    unsigned _nTotal;
    std::shared_ptr<spdlog::logger> _log;

  public:
    /**
     * constructor of class Closure
     * @param pointer to problem class
     */
    Closure( Settings* settings, QuadratureGrid* quadGrid );
    virtual ~Closure();
    static Closure* Create( Settings* settings, QuadratureGrid* quadGrid );

    const Matrix& GetPhiTildeWf() { return _phiTildeWf; }
};

#endif    // CLOSURE_H
