#ifndef UNIFORMSPARSEGRID_H
#define UNIFORMSPARSEGRID_H

#include "sparsegrid.h"

class UniformSparseGrid : public SparseGrid
{
  private:
    UniformSparseGrid() = delete;

    virtual unsigned order( unsigned level );
    virtual void createGrid();
    virtual void computeNodes( unsigned dim, unsigned level, unsigned nodeCount, unsigned** nodeOrder, double**& nodes );
    virtual void computeWeights(
        unsigned dim, unsigned level, unsigned nodeCount, unsigned unionCount, unsigned** nodeOrder, unsigned** unionIndex, double*& weights );
    // virtual void computeNodes1D(unsigned order );
    virtual std::vector<double> computeWeights1D( unsigned order );

  public:
    UniformSparseGrid( Settings* settings );
    UniformSparseGrid( Settings* settings, unsigned level );
    virtual ~UniformSparseGrid();
};

#endif    // UNIFORMSPARSEGRID_H
