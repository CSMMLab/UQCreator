#ifndef NOPLOT_H
#define NOPLOT_H

#include "plotengine.h"

class noPlot : public PlotEngine
{
    private:
  public:
    noPlot() = delete;
    noPlot( Problem* problem );
    virtual ~noPlot();
    static noPlot* Create( Problem* problem );
    virtual void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 ){}
    virtual void
    Plot1D( const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2 ){}
    virtual void Plot1D( const Vector& x1, const Vector& y1 ){}
    virtual void Plot1D( const Vector& x1,
                         const Vector& y1,
                         const Vector& x2,
                         const Vector& y2 ){}
};

#endif // NOPLOT_H
