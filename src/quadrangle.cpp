#include "quadrangle.h"

#include <cmath>

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
    Vector d1{_nodes[0]->coords[0] - _nodes[1]->coords[0], _nodes[0]->coords[1] - _nodes[1]->coords[1]};
    Vector d2{_nodes[1]->coords[0] - _nodes[2]->coords[0], _nodes[1]->coords[1] - _nodes[2]->coords[1]};
    Vector d3{_nodes[2]->coords[0] - _nodes[3]->coords[0], _nodes[2]->coords[1] - _nodes[3]->coords[1]};
    Vector d4{_nodes[3]->coords[0] - _nodes[0]->coords[0], _nodes[3]->coords[1] - _nodes[0]->coords[1]};

    double a = sqrt( pow( d1[0], 2 ) + pow( d1[1], 2 ) );
    double b = sqrt( pow( d2[0], 2 ) + pow( d2[1], 2 ) );
    double c = sqrt( pow( d3[0], 2 ) + pow( d3[1], 2 ) );
    double d = sqrt( pow( d4[0], 2 ) + pow( d4[1], 2 ) );
    double T = 0.5 * ( a + b + c + d );

    double alpha = acos( ( d4[0] * d1[0] + d4[1] * d1[1] ) / ( a * d ) );
    double beta  = acos( ( d2[0] * d3[0] + d2[1] * d3[1] ) / ( b * c ) );

    _area = sqrt( ( T - a ) * ( T - b ) * ( T - c ) * ( T - d ) - a * b * c * d * pow( cos( 0.5 * ( alpha + beta ) ), 2 ) );

    _center = Vector{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] + _nodes[3]->coords[0] ) / 4.0,
                     ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] + _nodes[3]->coords[1] ) / 4.0};
    _neighbors.resize( 4 );
    _neighborIDs.resize( 4 );
    SetupEdges();
}

Quadrangle::~Quadrangle() {}

void Quadrangle::SetupEdges() {
    Vector mid{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] + _nodes[3]->coords[0] ) / 4.0,
               ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] + _nodes[3]->coords[1] ) / 4.0};
    _minEdge = 1e10;
    _maxEdge = -1e10;
    for( unsigned i = 0; i < _N - 1; ++i ) {
        Node* A           = _nodes[i];
        Node* B           = _nodes[i + 1];
        double length     = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[1] - B->coords[1], 2 ) );
        Vector unitNormal = getOutwardNormal( A, B );
        unitNormal /= norm( unitNormal );
        Vector scaledNormal = length * unitNormal;

        _edges[i] = new Edge{A, B, length, unitNormal, scaledNormal};
        if( _minEdge > length ) _minEdge = length;
        if( _maxEdge < length ) _minEdge = length;
    }
    Node* A           = _nodes[_N - 1];
    Node* B           = _nodes[0];
    double length     = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[1] - B->coords[1], 2 ) );
    Vector unitNormal = getOutwardNormal( A, B );
    unitNormal /= norm( unitNormal );
    Vector scaledNormal = length * unitNormal;

    _edges[_N - 1] = new Edge{A, B, length, unitNormal, scaledNormal};
    if( _minEdge > length ) _minEdge = length;
    if( _maxEdge < length ) _minEdge = length;
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
