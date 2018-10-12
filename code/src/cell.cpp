#include "cell.h"

Cell::Cell( CELL_TYPE type, unsigned id, std::vector<Node*> nodes ) : _type( type ), _id( id ), _nodes( nodes ) {
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

Cell::~Cell() {}

Node* Cell::GetNode( unsigned i ) { return _nodes[i]; }

std::vector<Node*> Cell::GetNodes() { return _nodes; }

void Cell::AddNeighbor( const Cell* n ) {
    _neighbors.push_back( const_cast<Cell*>( n ) );
    _neighborIDs.resize( _neighbors.size() );
    _neighborIDs[_neighbors.size() - 1] = n->_id;
}

std::vector<Cell*> Cell::GetNeighbors() { return _neighbors; }

bool Cell::IsBoundaryCell() { return _isBoundaryCell; }

std::vector<Edge*> Cell::GetEdges() { return _edges; }

unsigned Cell::GetID() { return _id; }

double Cell::GetArea() { return _area; }

const Vector& Cell::GetCenter() { return _center; }

blaze::DynamicVector<unsigned> Cell::GetNeighborIDs() {}
