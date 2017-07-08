#ifndef MESH_H
#define MESH_H

#define MESH_1DPLAIN 100
#define MESH_MSH 101


typedef blaze::DynamicVector<double> vector;
typedef std::vector<blaze::DynamicVector<double> > meshData;

template <class T>
class Mesh
{
private:
    int _status;
    int _dimension;
    int _meshtype;
    int _numCells;
    T _mesh;
    T _spacing;
public:
    void load(std::string filename, int meshtype);
    void createGrid(double a, double b, int numCells);
    T* getGrid();
    T* getSpacing();
    vector
    Mesh();
};

#endif // MESH_H
