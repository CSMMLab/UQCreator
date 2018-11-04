#include "quadrangle.h"

Quadrangle::Quadrangle( unsigned id, std::vector<Node*> nodes ) : Cell( CELL_TYPE::QUADRILATERAL, id, nodes ) {
    unsigned boundaryNodeCtr = 0;
    for( const auto& n : _nodes ) {
        if( n->isBoundaryNode ) {
            boundaryNodeCtr++;
        }
    }
    if( boundaryNodeCtr == 2 ) {
        _isBoundaryCell = true;
    }
    else {
        _isBoundaryCell = false;
    }
    Vector d1{_nodes[0]->coords[0] - _nodes[2]->coords[0], _nodes[1]->coords[1] - _nodes[3]->coords[1]};
    Vector d2{_nodes[1]->coords[0] - _nodes[3]->coords[0], _nodes[1]->coords[1] - _nodes[3]->coords[1]};
    double angle = std::acos( dot( d1, d2 ) / ( norm( d1 ) * norm( d2 ) ) );

    _area = 0.5 * std::abs( dot( d1, d2 ) * std::sin( angle ) );

    _center = Vector{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] + _nodes[3]->coords[0] ) / 4.0,
                     ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] + _nodes[3]->coords[1] ) / 4.0};
    SetupEdges();
}

Quadrangle::~Quadrangle() {}

void Quadrangle::SetupEdges() {
    Vector mid{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] + _nodes[3]->coords[0] ) / 4.0,
               ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] + _nodes[3]->coords[1] ) / 4.0};
    for( unsigned i = 0; i < _N - 1; ++i ) {
        Node* A           = _nodes[i];
        Node* B           = _nodes[i + 1];
        double length     = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[1] - B->coords[1], 2 ) );
        Vector unitNormal = getOutwardNormal( A, B );
        unitNormal /= norm( unitNormal );
        Vector scaledNormal = length * unitNormal;

        _edges[i] = new Edge{A, B, length, unitNormal, scaledNormal};
    }
}

Vector Quadrangle::getOutwardNormal( Node* A, Node* B ) {
    double dx = A->coords[0] - B->coords[0];
    double dy = A->coords[1] - B->coords[1];
    Vector n{-dy, dx};
    Vector p{A->coords[0], A->coords[1]};
    Vector a{_center[0], _center[1]};
    if( dot( n, a - p ) > 0 ) {
        n *= -1;
    }
    n /= norm( n );
    return n;
}

unsigned Quadrangle::GetNodeNum() { return 4; }
