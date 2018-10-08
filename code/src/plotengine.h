#ifndef PLOTENGINE_H
#define PLOTENGINE_H

#include "blaze/math/DynamicVector.h"
#include <cpptoml.h>

#include "problem.h"

class PlotEngine
{
  private:
  protected:
    Problem* _problem;
    std::string _outputDir;
    std::vector<double> BlazeToStdVector( const Vector& v );

  public:
    PlotEngine() = delete;
    PlotEngine( Problem* problem );
    virtual ~PlotEngine();
    static PlotEngine* Create( Problem* problem );
    virtual void Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 ) = 0;
    virtual void
    Plot1D( const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2 ) = 0;
    virtual void Plot1D( const Vector& x1, const Vector& y1 )                                = 0;
    virtual void Plot1D( const Vector& x1,
                         const Vector& y1,
                         const Vector& x2,
                         const Vector& y2 )                                                                        = 0;
};

#endif    // PLOTENGINE_H
