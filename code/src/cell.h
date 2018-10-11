#ifndef ELEMENT_H
#define ELEMENT_H

#include <blaze/Math.h>
#include <iostream>
#include <vector>

#include "typedefs.h"

struct Node {
    unsigned id;
    bool isBoundaryNode;
    std::vector<double> coords;
};

struct Edge {
    Node* A;
    Node* B;
    double length;
    Vector unitNormal;
    Vector scaledNormal;
};

enum CELL_TYPE { TRIANGLE, QUADRILATERAL };

class Cell
{
  private:
    Cell();

  protected:
    CELL_TYPE _type;
    unsigned _id;
    unsigned _N;
    std::vector<Node*> _nodes;
    std::vector<Edge*> _edges;
    std::vector<Cell*> _neighbours;
    bool _isBoundaryCell;
    double _area;

    virtual void SetupEdges() = 0;

  public:
    Cell( CELL_TYPE type, unsigned pid, std::vector<Node*> pnodes );
    virtual ~Cell();

    Node* GetNode( unsigned i );
    std::vector<Node*> GetNodes();
    virtual unsigned GetNodeNum() = 0;
    void AddNeighbour( const Cell* n );
    std::vector<Cell*> GetNeighbours();
    bool IsBoundaryCell();
    std::vector<Edge*> GetEdges();
    unsigned GetID();
    double GetArea();
};

#endif    // ELEMENT_H
