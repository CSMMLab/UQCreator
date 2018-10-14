#ifndef MESH_H
#define MESH_H

#include <assert.h>
#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <iostream>
#include <string.h>

#include "cell.h"
#include "typedefs.h"

enum BoundaryType { NOSLIP, DIRICHLET, NEUMANN, PERIODIC, NONE };

struct BoundaryElement {
    unsigned type;
    std::vector<unsigned> nodes;
};

struct Boundary {
    std::string name;
    BoundaryType type;
    std::vector<BoundaryElement> elements;
};

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

    std::vector<blaze::DynamicVector<unsigned>> _neighborIDs;
    std::vector<BoundaryType> _boundaryType;

  public:
    static Mesh* Create( std::string inputFile );

    unsigned GetNumCells() const;
    unsigned GetDimension() const;
    std::vector<Cell*>& GetGrid();
    double GetArea( unsigned i ) const;
    Vector GetCenterPos( unsigned i ) const;

    unsigned GetNBoundaries() const;

    virtual std::vector<Cell*> GetNeighbors( unsigned i ) const;
    virtual blaze::DynamicVector<unsigned> GetNeighborIDs( unsigned i ) const;
    virtual Vector GetNormals( unsigned i, unsigned l ) const;
    virtual Vector GetUnitNormals( unsigned i, unsigned l ) const;

    BoundaryType GetBoundaryType( unsigned i );

    virtual Vector GetNodePositionsX() const = 0;

    virtual void Export( Matrix results ) const = 0;

    Mesh( unsigned dimension );
    virtual ~Mesh();
};

#endif    // MESH_H
