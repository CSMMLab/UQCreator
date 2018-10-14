#include "mesh.h"
#include "mesh1d.h"
#include "mesh2d.h"
// TODO: read Boundaries
Mesh::Mesh( unsigned dimension ) : _dimension( dimension ), _nBoundaries( 2 ) {}

Mesh::~Mesh() {}

Mesh* Mesh::Create( std::string inputFile ) {
    auto file     = cpptoml::parse_file( inputFile );
    auto settings = file->get_table( "mesh" );
    unsigned dim  = settings->get_as<unsigned>( "Dimension" ).value_or( 0 );
    if( dim == 1 ) {
        return new Mesh1D( inputFile );
    }
    else if( dim == 2 ) {
        return new Mesh2D( inputFile );
    }
    return nullptr;
}

unsigned Mesh::GetNumCells() const { return _numCells; }

unsigned Mesh::GetDimension() const { return _dimension; }

std::vector<Cell*>& Mesh::GetGrid() { return _cells; }

double Mesh::GetArea( unsigned i ) const { return _cells[i]->GetArea(); }

std::vector<Cell*> Mesh::GetNeighbors( unsigned i ) const { return _cells[i]->GetNeighbors(); }

Vector Mesh::GetNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->scaledNormal; }

Vector Mesh::GetUnitNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->unitNormal; }

Vector Mesh::GetCenterPos( unsigned i ) const { return _cells[i]->GetCenter(); }

unsigned Mesh::GetNBoundaries() const { return _nBoundaries; }

blaze::DynamicVector<unsigned> Mesh::GetNeighborIDs( unsigned i ) const { return _neighborIDs[i]; }

BoundaryType Mesh::GetBoundaryType( unsigned i ) { return _boundaryType[i]; }
