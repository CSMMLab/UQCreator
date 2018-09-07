#ifndef GNUPLOT_H
#define GNUPLOT_H

#include "plotengine.h"
#include <gnuplot-iostream.h>

class GnuplotLib : public PlotEngine
{
  private:
  public:
    GnuplotLib() = delete;
    GnuplotLib( Problem* problem );
    virtual ~GnuplotLib();
    void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 );
    void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2 );
    void Plot1D( const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1 );
    void Plot1D( const blaze::DynamicVector<double>& x1,
                 const blaze::DynamicVector<double>& y1,
                 const blaze::DynamicVector<double>& x2,
                 const blaze::DynamicVector<double>& y2 );
};

#endif    // GNUPLOT_H
