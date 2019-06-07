#ifndef TENSORIZEDQUADRATURE_H
#define TENSORIZEDQUADRATURE_H

#include <vector>

#include "legendre.h"
#include "quadraturegrid.h"
#include "settings.h"
#include "typedefs.h"

class TensorizedQuadrature : public QuadratureGrid
{
    const unsigned _nQuadPoints;
    std::vector<Polynomial*> _quad;
    Settings* _settings;
    unsigned _numDimXi;
    unsigned _nQTotal;
    double** _nodes;
    double* _weights;

  public:
    TensorizedQuadrature( Settings* settings );
    ~TensorizedQuadrature();
    std::vector<Vector> GetNodes();
    Vector GetWeights();
    unsigned GetNodeCount();
};

#endif    // TENSORIZEDQUADRATURE_H
