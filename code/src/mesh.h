#ifndef MESH_H
#define MESH_H

#include <string.h>
#include <assert.h>
#include <iostream>
#include <blaze/math/DynamicVector.h>
#include <cpptoml.h>

#define MESH_STATUS_UNLOADED 100
#define MESH_STATUS_LOADED 101
#define MESH_TYPE_1DPLAIN 110

class Mesh
{
private:
    int _status;
    int _dimension;
    int _meshType;
    int _numCells;
    blaze::DynamicVector<double> _mesh;
    blaze::DynamicVector<double> _spacing;

    Mesh(){}
public:
    void Load(std::string filename);
    void CreateGrid(double a, double b);

    int GetNumCells() const;
    int GetDimension() const;
    const blaze::DynamicVector<double>& GetGrid() const;
    const blaze::DynamicVector<double>& GetSpacing() const;

    Mesh(std::string inputFile);
    ~Mesh();
};

#endif // MESH_H
