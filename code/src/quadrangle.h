#ifndef QUADRANGLE_H
#define QUADRANGLE_H

#include "cell.h"

#include "assert.h"

class Quadrangle : public Cell
{

  private:
    Quadrangle();

    Vector getOutwardNormal( Node* A, Node* B );

  public:
    Quadrangle( unsigned id, std::vector<Node*> nodes );
    virtual ~Quadrangle();
    virtual void SetupEdges();
    virtual unsigned GetNodeNum();
};

#endif    // QUADRANGLE_H
