#ifndef LINE_H
#define LINE_H

#include "cell.h"

class Line : public Cell
{
private:
    Line();

public:
    Line( unsigned id, std::vector<Node*> nodes );
    virtual ~Line();
    virtual void SetupEdges();
    virtual unsigned GetNodeNum();
};

#endif // LINE_H
