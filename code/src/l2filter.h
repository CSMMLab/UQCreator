#ifndef L2FILTER_H
#define L2FILTER_H

#include "filter.h"

class L2Filter : public Filter
{
  private:
    L2Filter() = delete;
    double _lambda;
    Vector _filterFunction;

  public:
    L2Filter( Settings* settings );
    virtual ~L2Filter();
    void SetupFilter();
};

#endif    // L2FILTER_H
