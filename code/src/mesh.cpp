#include "mesh.h"
#include "mesh1d.h"
#include "mesh2d.h"

Mesh::Mesh( Settings* settings, unsigned dimension ) : _settings( settings ), _dimension( dimension ), _nBoundaries( 2 ) {
    _log        = spdlog::get( "event" );
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
    auto log     = spdlog::get( "event" );
    unsigned dim = settings->GetMeshDimension();
    if( dim == 1 ) {
        return new Mesh1D( settings );
    }
    else if( dim == 2 ) {
        return new Mesh2D( settings );
    }
    else {
        log->error( "[mesh] Unsupported mesh dimension: {0}", std::to_string( dim ) );
    }
    return nullptr;
}

unsigned Mesh::GetNumCells() const { return _numCells; }

unsigned Mesh::GetDimension() const { return _dimension; }

std::vector<Cell*>& Mesh::GetGrid() { return _cells; }

std::vector<Cell*> Mesh::GetGrid() const { return _cells; }

double Mesh::GetArea( unsigned i ) const { return _cells[i]->GetArea(); }

double Mesh::GetMinEdge( unsigned i ) const { return _cells[i]->GetMinEdge(); }

double Mesh::GetMaxEdge( unsigned i ) const { return _cells[i]->GetMaxEdge(); }

std::vector<Cell*> Mesh::GetNeighbors( unsigned i ) const { return _cells[i]->GetNeighbors(); }

Vector Mesh::GetNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->scaledNormal; }

Vector Mesh::GetUnitNormals( unsigned i, unsigned l ) const { return _cells[i]->GetEdges()[l]->unitNormal; }

Vector Mesh::GetCenterPos( unsigned i ) const { return _cells[i]->GetCenter(); }

unsigned Mesh::GetNBoundaries() const { return _nBoundaries; }

VectorU Mesh::GetNeighborIDs( unsigned i ) const { return _neighborIDs[i]; }

BoundaryType Mesh::GetBoundaryType( unsigned i ) const { return _boundaryType[i]; }

void Mesh::PlotInXi( const Matrix& u, unsigned state ) const {
    std::ofstream out( "../results/plotInXi" );
    unsigned nQ   = u.columns();
    unsigned nQ1D = unsigned( std::pow( nQ, 1.0 / _settings->GetNDimXi() ) );
    for( unsigned k = 0; k < nQ; ++k ) {
        if( k != 0 && k % nQ1D == 0 ) {
            out << std::endl;
        }
        out << u( k, state ) << " ";
    }
    out.close();
}

double Mesh::GetDomainArea() const {
    double area = 0;
    for( unsigned j = 0; j < _cells.size(); ++j ) {
        area += this->GetArea( j );
    }
    return area;
}
