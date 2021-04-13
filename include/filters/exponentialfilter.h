#ifndef EXPONENTIALFILTER_H
#define EXPONENTIALFILTER_H

#include "filter.h"

class ExponentialFilter : public Filter
{
  private:
    ExponentialFilter() = delete;
    double _c;    // machine precision constant

    double FilterFunction( double eta ) const;

  public:
    ExponentialFilter( Settings* settings );
    virtual ~ExponentialFilter();
    void SetupFilter();
};
#endif    // EXPONENTIALFILTER_H
