#ifndef ELEMENT_H
#define ELEMENT_H

#include <assert.h>
#include <iostream>
#include <vector>

#include "typedefs.h"

enum BoundaryType { NOSLIP, DIRICHLET, NEUMANN, PERIODIC, NONE, SWWALL };

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

enum CELL_TYPE { LINE, TRIANGLE, QUADRILATERAL };

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
    std::vector<Cell*> _neighbors;
    VectorU _neighborIDs;
    bool _isBoundaryCell;
    BoundaryType _boundaryType;
    double _area;
    Vector _center;
    Vector _boundaryNormal;
    double _minEdge;    // saves edge with minimal length for CFL computation
    double _maxEdge;    // saves edge with minimal length for CFL computation

    virtual void SetupEdges() = 0;

  public:
    Cell( CELL_TYPE type, unsigned pid, std::vector<Node*> pnodes );
    virtual ~Cell();

    Node* GetNode( unsigned i );
    std::vector<Node*> GetNodes();
    virtual unsigned GetNodeNum() = 0;
    void AddNeighbor( const Cell* n, unsigned k );
    void AddNeighborId( unsigned n, unsigned k );
    std::vector<Cell*> GetNeighbors();
    VectorU GetNeighborIDs();
    bool IsBoundaryCell();
    std::vector<Edge*> GetEdges();
    unsigned GetNPoints() const;
    unsigned GetID();
    double GetArea();
    double GetMinEdge() const;
    double GetMaxEdge() const;
    const Vector& GetCenter();
    void SetBoundaryType( BoundaryType type );
    BoundaryType GetBoundaryType() const;
    void UpdateBoundaryNormal();
    Vector GetBoundaryUnitNormal();
    Vector GetUnitNormal( unsigned i );
    Vector GetNormal( unsigned i );
    void SetDefaultCellId( unsigned id );
};

#endif    // ELEMENT_H
