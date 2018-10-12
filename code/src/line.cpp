#include "line.h"

Line::Line( unsigned id, std::vector<Node*> nodes ) : Cell( CELL_TYPE::LINE, id, nodes ) {
    unsigned boundaryNodeCtr = 0;
    for( const auto& n : _nodes ) {
        if( n->isBoundaryNode ) {
            boundaryNodeCtr++;
        }
    }
    if( boundaryNodeCtr == 1 ) {
        _isBoundaryCell = true;
    }
    else {
        _isBoundaryCell = false;
    }
    _area   = std::fabs( _nodes[1]->coords[0] - _nodes[0]->coords[0] );
    _center = Vector{( _nodes[0]->coords[0] + _nodes[0]->coords[0] ) / 2};

    SetupEdges();
}

Line::~Line() {}

void Line::SetupEdges() {
    Node *A, *B;
    if( _nodes[0]->coords[0] < _nodes[1]->coords[0] ) {
        A = _nodes[0];
        B = _nodes[1];
    }
    else {
        B = _nodes[0];
        A = _nodes[1];
    }
    double length = 1;
    Vector unitNormal{1};
    Vector scaledNormal = length * unitNormal;

    _edges[0] = new Edge{A, A, 1, -unitNormal, -scaledNormal};
    _edges[1] = new Edge{B, B, 1, unitNormal, scaledNormal};
}

unsigned Line::GetNodeNum() { return 2; }
