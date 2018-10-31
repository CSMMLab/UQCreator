#include "mesh.h"
#include "mesh1d.h"
#include "mesh2d.h"

Mesh::Mesh( Settings* settings, unsigned dimension ) : _settings( settings ), _dimension( dimension ), _nBoundaries( 2 ) {
    _outputFile = _settings->GetOutputFile();
}

Mesh::~Mesh() {
    for( auto& c : _cells ) {
        delete c;
    }
    for( auto& n : _nodes ) {
        delete n;
    }
}

Mesh* Mesh::Create( Settings* settings ) {
    unsigned dim = settings->GetMeshDimension();
    if( dim == 1 ) {
        return new Mesh1D( settings );
    }
    else if( dim == 2 ) {
        return new Mesh2D( settings );
    }
    else {
        std::cerr << "[Mesh] Unsupported mesh dimension: " + std::to_string( dim ) << std::endl;
    }
    return nullptr;
}

unsigned Mesh::GetNumCells() const { return _numCells; }

unsigned Mesh::GetDimension() const { return _dimension; }

std::vector<Cell*>& Mesh::GetGrid() { return _cells; }

std::vector<Cell*> Mesh::GetGrid() const { return _cells; }

double Mesh::GetArea( unsigned i ) const { return _cells[i]->GetArea(); }

std::vector<Cell*> Mesh::GetNeighbors( unsigned i ) const { return _cells[i]->GetNeighbors(); }

Vector Mesh::GetNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->scaledNormal; }

Vector Mesh::GetUnitNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->unitNormal; }

Vector Mesh::GetCenterPos( unsigned i ) const { return _cells[i]->GetCenter(); }

unsigned Mesh::GetNBoundaries() const { return _nBoundaries; }

VectorU Mesh::GetNeighborIDs( unsigned i ) const { return _neighborIDs[i]; }

BoundaryType Mesh::GetBoundaryType( unsigned i ) const { return _boundaryType[i]; }
