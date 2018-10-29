#ifndef MESH_H
#define MESH_H

#include <assert.h>
#include <cpptoml.h>
#include <iostream>
#include <string.h>

#include "cell.h"
#include "settings.h"
#include "typedefs.h"

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
    Mesh() = delete;

  protected:
    Settings* _settings;
    unsigned _dimension;
    unsigned _numCells;
    unsigned _nBoundaries;
    std::vector<Cell*> _cells;
    std::vector<Node*> _nodes;
    std::string _outputFile;

    std::vector<VectorU> _neighborIDs;
    std::vector<BoundaryType> _boundaryType;

  public:
    static Mesh* Create( Settings* settings );

    unsigned GetNumCells() const;
    unsigned GetDimension() const;
    std::vector<Cell*>& GetGrid();
    std::vector<Cell*> GetGrid() const;
    double GetArea( unsigned i ) const;
    Vector GetCenterPos( unsigned i ) const;

    unsigned GetNBoundaries() const;

    virtual std::vector<Cell*> GetNeighbors( unsigned i ) const;
    virtual VectorU GetNeighborIDs( unsigned i ) const;
    virtual Vector GetNormals( unsigned i, unsigned l ) const;
    virtual Vector GetUnitNormals( unsigned i, unsigned l ) const;

    BoundaryType GetBoundaryType( unsigned i ) const;

    virtual Vector GetNodePositionsX() const = 0;

    virtual void Export( Matrix results ) const = 0;

    Mesh( Settings* settings, unsigned dimension );
    virtual ~Mesh();
};

#endif    // MESH_H
