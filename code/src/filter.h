#ifndef FILTER_H
#define FILTER_H

#include <spdlog/spdlog.h>
#include <vector>

#include "settings.h"
#include "typedefs.h"

class Filter
{
  protected:
    Settings* _settings;
    Filter() = delete;
    double _lambda;    // filter strength
    Vector _filterFunction;
    unsigned _filterOrder;    // order of the filter
    unsigned _gamma;
    std::shared_ptr<spdlog::logger> _log;

    double FilterFunction( double eta ) const;

  public:
    Filter( Settings* settings );
    virtual ~Filter();
    virtual void SetupFilter();
    virtual void FilterMoments( Tensor& u ) const;
    virtual void FilterMoments( Tensor& v, const Tensor& u ) const;
    virtual void FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const;
    static Filter* Create( Settings* settings );
};

#endif    // FILTER_H
