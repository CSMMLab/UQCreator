#ifndef MESH2D_H
#define MESH2D_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include "triangle.h"

using vtkPointsSP                    = vtkSmartPointer<vtkPoints>;
using vtkUnstructuredGridSP          = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkTriangleSP                  = vtkSmartPointer<vtkTriangle>;
using vtkCellArraySP                 = vtkSmartPointer<vtkCellArray>;
using vtkXMLUnstructuredGridWriterSP = vtkSmartPointer<vtkXMLUnstructuredGridWriter>;

struct BoundaryElement {
    unsigned type;
    std::vector<unsigned> nodes;
};

struct Boundary {
    std::string name;
    std::vector<BoundaryElement> elements;
};

class Mesh2D
{
  private:
    unsigned _dim;
    std::vector<Cell*> _elements;
    std::vector<Node*> _nodes;
    std::vector<Boundary> _boundaries;

    unsigned GetTrailingPosNumber( std::string s );
    unsigned BinarySearch( unsigned id, unsigned lBound, unsigned rBound );
    Node* FindNodeByID( unsigned id );
    void DetermineNeighbours();

  public:
    Mesh2D();
    ~Mesh2D();
    void LoadSU2MeshFromFile( std::string meshfile );
    void ExportToVTK( std::string vtkfile );
    double GetArea( unsigned i );
};

#endif    // MESH2D_H
