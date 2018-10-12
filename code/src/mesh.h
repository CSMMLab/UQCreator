#ifndef MESH_H
#define MESH_H

#include <assert.h>
#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <iostream>
#include <string.h>

#include "cell.h"
#include "typedefs.h"

class Mesh
{
  private:
    Mesh() {}

  protected:
    unsigned _dimension;
    unsigned _numCells;
    unsigned _nBoundaries;
    std::vector<Cell*> _cells;
    std::vector<Node*> _nodes;
    std::string _outputFile;

    std::vector<blaze::DynamicVector<unsigned>> _neighbors;

  public:
    static Mesh* Create( std::string inputFile );

    unsigned GetNumCells() const;
    unsigned GetDimension() const;
    std::vector<Cell*>& GetGrid();
    double GetArea( unsigned i ) const;
    Vector GetCenterPos( unsigned i ) const;

    unsigned GetNBoundaries() const;

    virtual std::vector<Cell*> GetNeighbors( unsigned i ) const;
    virtual blaze::DynamicVector<unsigned> GetNeighborsIndex( unsigned i ) const;
    virtual Vector GetNormals( unsigned i, unsigned l ) const;
    virtual Vector GetUnitNormals( unsigned i, unsigned l ) const;

    virtual Vector GetNodePositionsX() const {};

    virtual void Export() const = 0;

    Mesh( unsigned dimension );
    virtual ~Mesh();
};

#endif    // MESH_H
