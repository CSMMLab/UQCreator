#ifndef HOULIFILTER_H
#define HOULIFILTER_H

#include "filter.h"

class HouLiFilter : public Filter
{
  private:
    HouLiFilter() = delete;

    double _eps;
    unsigned _gamma;

    double FilterFunction( double eta ) const;

  public:
    HouLiFilter( Settings* settings );
    void SetupFilter();
    virtual ~HouLiFilter();
};

#endif    // HOULIFILTER_H
