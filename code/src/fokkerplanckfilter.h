#ifndef FOKKERPLANCKFILTER_H
#define FOKKERPLANCKFILTER_H

#include "filter.h"

class FokkerPlanckFilter : public Filter
{
  private:
    FokkerPlanckFilter() = delete;
    double _lambda;
    Vector _filterFunction;

  public:
    FokkerPlanckFilter( Settings* settings );
    virtual ~FokkerPlanckFilter();
    void SetupFilter();
};

#endif // FOKKERPLANCKFILTER_H
