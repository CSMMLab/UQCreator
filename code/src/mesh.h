#ifndef MESH_H
#define MESH_H

#include <assert.h>
#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>
#include <iostream>
#include <string.h>

#include "typedefs.h"

#define MESH_STATUS_UNLOADED 100
#define MESH_STATUS_LOADED 101
#define MESH_TYPE_1DPLAIN 110

class Mesh
{
  private:
    int _status;
    int _meshType;
    unsigned _dimension;
    unsigned _numCells;
    Vector _mesh;
    Vector _spacing;
    std::vector<blaze::DynamicVector<unsigned>> _neighbors;
    std::vector<std::vector<Vector>> _n;
    std::vector<std::vector<Vector>> _nUnit;

    Mesh() {}

  public:
    void Load( std::string filename );
    void CreateGrid( double a, double b );

    unsigned GetNumCells() const;
    unsigned GetDimension() const;
    const Vector& GetGrid() const;
    const Vector& GetSpacing() const;
    double GetArea( unsigned j ) const;
    blaze::DynamicVector<unsigned> GetNeighbors( unsigned j ) const;
    Vector GetNormals( unsigned j, unsigned l ) const;
    Vector GetUnitNormals( unsigned j, unsigned l ) const;

    Mesh( std::string inputFile );
    ~Mesh();
};

#endif    // MESH_H
