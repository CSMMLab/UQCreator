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
    void load(std::string filename);
    void createGrid(double a, double b);

    int getNumCells() const;
    int getDimension() const;
    meshData& getGrid();
    meshData& getSpacing();

    Mesh(std::string inputFile);
};

#endif // MESH_H
