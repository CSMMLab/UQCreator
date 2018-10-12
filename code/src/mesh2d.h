#ifndef MESH2D_H
#define MESH2D_H

#include <algorithm>
#include <assert.h>
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

#include "mesh.h"
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

enum BoundaryType { NOSLIP, DIRICHLET, NEUMANN, PERIODIC };

class Mesh2D : public Mesh
{
  private:
    std::vector<Boundary> _boundaries;
    std::string _SU2MeshFile;

    unsigned GetTrailingNumber( std::string s );
    unsigned BinarySearch( unsigned id, unsigned lBound, unsigned rBound );
    Node* FindNodeByID( unsigned id );
    void DetermineNeighbors();
    void LoadSU2MeshFromFile( std::string meshfile );
    void ExportToVTK( std::string vtkfile ) const;

    Mesh2D();

  public:
    Mesh2D( std::string inputFile );
    virtual ~Mesh2D();

    virtual void Export() const;
};

#endif    // MESH2D_H
