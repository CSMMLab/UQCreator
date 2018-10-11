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
    }
    else {
        _isBoundaryCell = false;
    }
    _area = std::abs( ( _nodes[0]->coords[0] * ( _nodes[1]->coords[1] - _nodes[2]->coords[1] ) +
                        _nodes[1]->coords[0] * ( _nodes[2]->coords[1] - _nodes[0]->coords[1] ) +
                        _nodes[2]->coords[0] * ( _nodes[0]->coords[1] - _nodes[1]->coords[1] ) ) /
                      2 );
}

Triangle::~Triangle() {}

void Triangle::SetupEdges() {
    Vector mid{( _nodes[0]->coords[0] + _nodes[1]->coords[0] + _nodes[2]->coords[0] ) / 3,
               ( _nodes[0]->coords[1] + _nodes[1]->coords[1] + _nodes[2]->coords[1] ) / 3};
    for( unsigned i = 0; i < _N - 1; ++i ) {
        Node* A       = _nodes[i];
        Node* B       = _nodes[i + 1];
        double length = std::sqrt( std::pow( A->coords[0] - B->coords[0], 2 ) + std::pow( A->coords[0] - B->coords[0], 2 ) );
        double dx     = A->coords[0] - B->coords[0];
        double dy     = A->coords[1] - B->coords[1];
        Vector unitNormal{dy, -dx};
        if( unitNormal == midNormal( mid, A, B ) ) {
            unitNormal[0] = -dy;
            unitNormal[1] = dx;
        }
        assert( unitNormal == midNormal( mid, A, B ) );
        unitNormal /= blaze::norm( unitNormal );
        Vector scaledNormal = length * unitNormal;
        _edges.push_back( new Edge{A, B, length, unitNormal, scaledNormal} );
    }
}

Vector Triangle::midNormal( Vector mid, Node* A, Node* B ) {
    Vector midLine{( A->coords[0] + B->coords[0] ) / 2, ( A->coords[1] + B->coords[1] ) / 2};
    Vector res{midLine[0] - mid[0], midLine[1] - mid[1]};
    res /= norm( res );
    return res;
}

unsigned Triangle::GetNodeNum() { return 3; }
