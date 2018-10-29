#include "cell.h"

Cell::Cell( CELL_TYPE type, unsigned id, std::vector<Node*> nodes ) : _type( type ), _id( id ), _nodes( nodes ), _boundaryType( BoundaryType::NONE ) {
    if( _type == CELL_TYPE::LINE ) {
        _N = 2;
    }
    else if( _type == CELL_TYPE::TRIANGLE ) {
        _N = 3;
    }
    else if( _type == CELL_TYPE::QUADRILATERAL ) {
        _N = 4;
    }
    _edges.resize( _N );
}

Cell::~Cell() {
    for( auto& e : _edges ) {
        delete e;
    }
}

Node* Cell::GetNode( unsigned i ) { return _nodes[i]; }

std::vector<Node*> Cell::GetNodes() { return _nodes; }

void Cell::AddNeighbor( const Cell* n, unsigned k ) {
    _neighbors[k]   = const_cast<Cell*>( n );
    _neighborIDs[k] = n->_id;
}

void Cell::AddNeighborId( unsigned n, unsigned k ) { _neighborIDs[k] = n; }

std::vector<Cell*> Cell::GetNeighbors() { return _neighbors; }

bool Cell::IsBoundaryCell() { return _isBoundaryCell; }

std::vector<Edge*> Cell::GetEdges() { return _edges; }

unsigned Cell::GetID() { return _id; }

double Cell::GetArea() { return _area; }

const Vector& Cell::GetCenter() { return _center; }

VectorU Cell::GetNeighborIDs() { return _neighborIDs; }

void Cell::SetBoundaryType( BoundaryType type ) { _boundaryType = type; }

BoundaryType Cell::GetBoundaryType() const { return _boundaryType; }

Vector Cell::GetBoundaryUnitNormal() { return _boundaryNormal; }

void Cell::UpdateBoundaryNormal() {
    assert( this->IsBoundaryCell() );
    auto elemPtr    = std::max_element( _neighborIDs.begin(), _neighborIDs.end() );
    unsigned id     = static_cast<unsigned>( std::distance( _neighborIDs.begin(), elemPtr ) );
    _boundaryNormal = _edges[id]->unitNormal;
}

Vector Cell::GetUnitNormal( unsigned i ) { return _edges[i]->unitNormal; }

Vector Cell::GetNormal( unsigned i ) { return _edges[i]->scaledNormal; }
