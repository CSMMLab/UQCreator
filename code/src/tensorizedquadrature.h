#ifndef TENSORIZEDQUADRATURE_H
#define TENSORIZEDQUADRATURE_H

#include <vector>

#include "legendre.h"
#include "quadraturegrid.h"
#include "settings.h"
#include "typedefs.h"

class TensorizedQuadrature : public QuadratureGrid
{
    Settings* _settings;
    const unsigned _nQuadPoints;
    std::vector<Polynomial*> _quad;    
    unsigned _numDimXi;
    unsigned _nQTotal;
    double** _nodes;
    double* _weights;
    void CreateGrid();

  public:
    TensorizedQuadrature( Settings* settings );
    TensorizedQuadrature( Settings* settings, unsigned numDimXi, unsigned nQuadPoints );
    ~TensorizedQuadrature();
    std::vector<Vector> GetNodes();
    Vector GetWeights();
    unsigned GetNodeCount();
};

#endif    // TENSORIZEDQUADRATURE_H
