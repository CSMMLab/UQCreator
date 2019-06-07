#ifndef QUADRATUREGRID_H
#define QUADRATUREGRID_H

#include "settings.h"
#include "typedefs.h"
#include <spdlog/spdlog.h>

class QuadratureGrid
{
private:
    Settings* _settings;
public:
    QuadratureGrid();
    virtual ~QuadratureGrid(){}
    static QuadratureGrid* Create( Settings* settings );

    virtual unsigned GetNodeCount() = 0;
    virtual std::vector<Vector> GetNodes() = 0;
    virtual Vector GetWeights() = 0;
};

#endif // QUADRATUREGRID_H
