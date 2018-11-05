#ifndef Mesh1D_H
#define Mesh1D_H

#include <assert.h>
#include <cpptoml.h>
#include <iostream>
#include <string.h>

#include "line.h"
#include "mesh.h"
#include "typedefs.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

using vtkPointsSP                    = vtkSmartPointer<vtkPoints>;
using vtkUnstructuredGridSP          = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkCellArraySP                 = vtkSmartPointer<vtkCellArray>;
using vtkDoubleArraySP               = vtkSmartPointer<vtkDoubleArray>;
using vtkXMLUnstructuredGridWriterSP = vtkSmartPointer<vtkXMLUnstructuredGridWriter>;
using vtkCellDataToPointDataSP       = vtkSmartPointer<vtkCellDataToPointData>;
using vtkQuadSP                      = vtkSmartPointer<vtkQuad>;

class Mesh1D : public Mesh
{
  private:
    int _status;
    int _MeshType;

    Mesh1D();
    void CreateGrid( double a, double b );

  public:
    virtual void Export( const Matrix& results ) const;

    Mesh1D( Settings* settings );
    virtual Vector GetNodePositionsX() const;
    virtual ~Mesh1D();
};

#endif    // Mesh1D_H
