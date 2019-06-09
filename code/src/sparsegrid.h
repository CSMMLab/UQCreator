#ifndef SPARSEGRID_H
#define SPARSEGRID_H

#include <assert.h>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#include "mathtools.h"
#include "quadraturegrid.h"
#include "settings.h"
#include "typedefs.h"

class SparseGrid : public QuadratureGrid
{
  protected:
    std::string _nodeType;
    unsigned _dim;
    unsigned _level;
    double** _nodes;
    double* _weights;
    unsigned _nodeCount;
    unsigned _unionCount;
    unsigned** _nodeOrder;
    unsigned** _unionIndex;
    std::vector<std::vector<std::vector<unsigned>>> _gridIndices;

    void computeUnionIndex( unsigned dim, unsigned level, unsigned& nodeCount, unsigned& _unionCount, unsigned**& unionIndex );
    void computeNodeOrder( unsigned dim,
                           unsigned level,
                           unsigned unionCount,
                           unsigned** unionIndex,
                           unsigned nodeCount,
                           std::vector<std::vector<std::vector<unsigned>>>& gridIndices,
                           unsigned**& nodeOrder );
    template <class T> void colexOrdering( std::vector<T> base, std::vector<T>& mat, bool& run );
    std::vector<double> prod( unsigned index, unsigned order, unsigned orderND, std::vector<double> weights1D, std::vector<double> weightsND );
    int modu( unsigned i );
    int choose( int n, int k );

    SparseGrid() = delete;

    virtual unsigned order( unsigned level )                                                                             = 0;
    virtual void createGrid()                                                                                            = 0;
    virtual void computeNodes( unsigned dim, unsigned level, unsigned nodeCount, unsigned** nodeOrder, double**& nodes ) = 0;
    virtual void computeWeights(
        unsigned dim, unsigned level, unsigned nodeCount, unsigned unionCount, unsigned** nodeOrder, unsigned** unionIndex, double*& weights ) = 0;
    // virtual void computeNodes1D(unsigned order ) = 0;
    virtual std::vector<double> computeWeights1D( unsigned order ) = 0;

  public:
    SparseGrid( unsigned dim, unsigned level );
    virtual ~SparseGrid();
    virtual unsigned GetNodeCount();
    virtual std::vector<Vector> GetNodes();
    virtual Vector GetWeights();
    std::vector<Vector> GetNodeOrder();
    std::vector<Vector> GetUnionIndex();
    unsigned GetUnionCount() { return _unionCount; }
    // static SparseGrid* Create( std::string type, unsigned dim, unsigned level );
};

#endif    // SPARSEGRID_H
