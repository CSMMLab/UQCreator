#ifndef ERFCLOGFILTER_H
#define ERFCLOGFILTER_H

#include "filter.h"

class ErfcLogFilter : public Filter
{
  private:
    ErfcLogFilter() = delete;

    double FilterFunction( double eta ) const;

  public:
    ErfcLogFilter( Settings* settings );
    virtual ~ErfcLogFilter();
    void SetupFilter();
};

#endif // ERFCLOGFILTER_H
