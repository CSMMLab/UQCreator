#ifndef SPLINEFILTER_H
#define SPLINEFILTER_H

#include "filter.h"

class SplineFilter : public Filter
{
  private:
    SplineFilter() = delete;
    double _eta;    // variance correction

    double FilterFunction( double eta ) const;

  public:
    SplineFilter( Settings* settings );
    virtual ~SplineFilter();
    void SetupFilter();
};

#endif    // SPLINEFILTER_H
