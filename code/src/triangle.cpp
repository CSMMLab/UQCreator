#include "triangle.h"

Triangle::Triangle( unsigned id, std::vector<Node*> nodes ) : Cell( CELL_TYPE::TRIANGLE, id, nodes ) {
    unsigned boundaryNodeCtr = 0;
    for( const auto& n : _nodes ) {
        if( n->isBoundaryNode ) {
            boundaryNodeCtr++;
        }
    }
    if( boundaryNodeCtr == 2 ) {
        _isBoundaryCell = true;
        // this->AddNeighbor(_numCells)
    }
    else {
        _isBoundaryCell = false;
    }
    _area   = std::abs( ( _nodes[0]->coords[0] * ( _nodes[1]->coords[1] - _nodes[2]->coords[1] ) +
                        _nodes[1]->coords[0] * ( _nodes[2]->coords[1] - _nodes[0]->coords[1] ) +
                        _nodes[2]->coords[0] * ( _nodes[0]->coords[1] - _nodes[1]->coords[1] ) ) /
                      2 );
    _center = Vector{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] ) / 3,
                     ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] ) / 3};
    SetupEdges();
    _neighbors.resize( 3 );
    _neighborIDs.resize( 3 );
}

Triangle::~Triangle() {}

void Triangle::SetupEdges() {
    Vector mid{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] ) / 3,
               ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] ) / 3};
    _minEdge = 1e10;
    for( unsigned i = 0; i < _N - 1; ++i ) {
        Node* A           = _nodes[i];
        Node* B           = _nodes[i + 1];
        double length     = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[1] - B->coords[1], 2 ) );
        Vector unitNormal = getOutwardNormal( A, B );
        unitNormal /= norm( unitNormal );
        Vector scaledNormal = length * unitNormal;

        _edges[i] = new Edge{A, B, length, unitNormal, scaledNormal};
        if( _minEdge > length ) _minEdge = length;
    }
    Node* A           = _nodes[_N - 1];
    Node* B           = _nodes[0];
    double length     = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[1] - B->coords[1], 2 ) );
    Vector unitNormal = getOutwardNormal( A, B );
    unitNormal /= norm( unitNormal );
    Vector scaledNormal = length * unitNormal;

    _edges[_N - 1] = new Edge{A, B, length, unitNormal, scaledNormal};
    if( _minEdge > length ) _minEdge = length;
}

Vector Triangle::getOutwardNormal( Node* A, Node* B ) {
    double dx = A->coords[0] - B->coords[0];
    double dy = A->coords[1] - B->coords[1];
    Vector n{-dy, dx};
    Vector p{A->coords[0], A->coords[1]};
    Vector a{_center[0], _center[1]};
    if( dot( n, a - p ) > 0 ) {
        n *= -1.0;
    }
    n /= norm( n );
    return n;
}

unsigned Triangle::GetNodeNum() { return 3; }
