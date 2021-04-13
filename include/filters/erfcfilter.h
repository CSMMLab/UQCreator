#ifndef ERFCFILTER_H
#define ERFCFILTER_H

#include "filter.h"

class ErfcFilter : public Filter
{
  private:
    ErfcFilter() = delete;

    double FilterFunction( double eta ) const;

  public:
    ErfcFilter( Settings* settings );
    virtual ~ErfcFilter();
    void SetupFilter();
};

#endif // ERFCFILTER_H
