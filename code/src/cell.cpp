#include "cell.h"

Cell::Cell( CELL_TYPE type, unsigned id, std::vector<Node*> nodes ) : _type( type ), _id( id ), _nodes( nodes ) {
    if( _type == CELL_TYPE::TRIANGLE ) {
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

void Cell::AddNeighbour( const Cell* n ) { _neighbours.push_back( const_cast<Cell*>( n ) ); }

std::vector<Cell*> Cell::GetNeighbours() { return _neighbours; }

bool Cell::IsBoundaryCell() { return _isBoundaryCell; }

std::vector<Edge*> Cell::GetEdges() { return _edges; }

unsigned Cell::GetID() { return _id; }
