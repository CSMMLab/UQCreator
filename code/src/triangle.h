#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "assert.h"
#include "cell.h"

class Triangle : public Cell
{
  private:
    Triangle();

    Vector midNormal( Node* A, Node* B );

  public:
    Triangle( unsigned id, std::vector<Node*> nodes );
    virtual ~Triangle();
    virtual void SetupEdges();
    virtual unsigned GetNodeNum();
};

#endif    // TRIANGLE_H
