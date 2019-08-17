#ifndef TENSORIZEDCC_H
#define TENSORIZEDCC_H

#include <vector>

#include "legendre.h"
#include "quadraturegrid.h"
#include "settings.h"
#include "typedefs.h"

class TensorizedCC : public QuadratureGrid
{
    Settings* _settings;
    const unsigned _nQuadPoints;
    std::vector<Polynomial*> _quad;
    unsigned _numDimXi;
    unsigned _nQTotal;
    double** _nodes;
    double* _weights;
    std::vector<std::vector<unsigned>> _counter;
    void CreateGrid();
    Vector computeNodes1D( unsigned level );
    Vector computeWeights1D( unsigned level );
    void DetermineCounter( unsigned level );

  public:
    TensorizedCC( Settings* settings );
    TensorizedCC( Settings* settings, unsigned level );
    ~TensorizedCC();
    std::vector<Vector> GetNodes();
    Vector GetWeights();
    unsigned GetNodeCount();
};

#endif    // TENSORIZEDCC_H
