#ifndef MESH2D_H
#define MESH2D_H

#include <algorithm>
#include <assert.h>
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
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
using vtkDoubleArraySP               = vtkSmartPointer<vtkDoubleArray>;
using vtkXMLUnstructuredGridWriterSP = vtkSmartPointer<vtkXMLUnstructuredGridWriter>;
using vtkCellDataToPointDataSP       = vtkSmartPointer<vtkCellDataToPointData>;

class Mesh2D : public Mesh
{
  private:
    enum MeshFormat { VTK };
    std::vector<std::pair<std::string, BoundaryType>> _BCs;
    MeshFormat _format;
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
    Mesh2D( Settings* settings );
    virtual ~Mesh2D();

    virtual Vector GetNodePositionsX() const;
    virtual void Export( Matrix results ) const;
};

#endif    // MESH2D_H
