#ifndef PLOTENGINE_H
#define PLOTENGINE_H

#include "blaze/math/DynamicVector.h"
#include <cpptoml.h>

class PlotEngine
{
private:

protected:
    std::vector<double> BlazeToStdVector(const blaze::DynamicVector<double>& v);

public:
    PlotEngine();
    virtual ~PlotEngine();
    static PlotEngine* Create(std::string inputFile);
    virtual void Plot1D(const std::vector<double>& x1, const std::vector<double>& y1) = 0;
    virtual void Plot1D(const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2) = 0;
    virtual void Plot1D(const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1) = 0;
    virtual void Plot1D(const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1, const blaze::DynamicVector<double>& x2, const blaze::DynamicVector<double>& y2) = 0;
};

#endif // PLOTENGINE_H
