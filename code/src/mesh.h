#ifndef MESH_H
#define MESH_H

#include <string.h>
#include <blaze/Blaze.h>
#include <cpptoml.h>

typedef blaze::DynamicVector<double> vector;
typedef std::vector<vector> meshData;

#define MESH_STATUS_UNLOADED 100
#define MESH_STATUS_LOADED 101
#define MESH_TYPE_1DPLAIN 110

class Mesh
{
private:
    int _status;
    int _dimension;
    int _meshtype;
    int _numCells;
    meshData _mesh;
    meshData _spacing;

    Mesh(){}
public:
    void Load(std::string filename);
    void CreateGrid(double a, double b);

    int GetNumCells() const;
    int GetDimension() const;
    meshData& GetGrid();
    meshData& GetSpacing();

    Mesh(std::string inputFile);
};

#endif // MESH_H
