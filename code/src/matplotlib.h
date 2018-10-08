#ifndef MATPLOTLIB_H
#define MATPLOTLIB_H

#include <chrono>
#include <matplotlibcpp.h>

#include "plotengine.h"

class Matplotlib : public PlotEngine
{
  private:
    const std::string _fileExt = ".pdf";

  public:
    Matplotlib() = delete;
    Matplotlib( Problem* problem );
    virtual ~Matplotlib();
    void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 );
    void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2 );
    void Plot1D( const Vector& x1, const Vector& y1 );
    void Plot1D( const Vector& x1,
                 const Vector& y1,
                 const Vector& x2,
                 const Vector& y2 );
};

#endif    // MATPLOTLIB_H
