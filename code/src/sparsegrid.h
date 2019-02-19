#ifndef SPARSE_GRID_H
#define SPARSE_GRID_H

#include <assert.h>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

class SparseGrid
{
  private:
    std::string _nodeType;
    unsigned _dim;
    unsigned _level;
    unsigned _Rl;
    long unsigned _nodeCount;
    long unsigned _unionCount;

    unsigned** _unionIndex;
    unsigned** _nodeOrder;
    std::vector<std::vector<std::vector<unsigned>>> _gridIndices;

    double** _nodes;
    double* _weights;

    void computeUnionIndex();
    void computeNodeOrder();
    void createCCgrid();
    void computeCCNodes();
    void computeCCWeights();
    std::vector<double> CCweights1D( unsigned order );
    std::vector<double>
    CCweightsND( unsigned dimIndex, unsigned order, unsigned orderND, std::vector<double> weights1D, std::vector<double> weightsND );
    void colexOrdering( std::vector<unsigned> base, std::vector<unsigned>& mat, bool& run );
    template <class T> int isEven( T i );
    template <class T> double choose( T n, T k );

    SparseGrid() = delete;

  public:
    SparseGrid( std::string type, unsigned pDim, unsigned pLevel );
    ~SparseGrid();
    std::vector<std::pair<std::vector<double>, double>> GetGrid() const;
    std::vector<std::vector<double>> GetNodes() const;
    std::vector<double> GetWeights() const;
    unsigned long GetNodeCount() const;
};

#endif
