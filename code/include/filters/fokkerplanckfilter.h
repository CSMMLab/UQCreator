#ifndef FOKKERPLANCKFILTER_H
#define FOKKERPLANCKFILTER_H

#include "filter.h"

class FokkerPlanckFilter : public Filter
{
  private:
    FokkerPlanckFilter() = delete;

  public:
    FokkerPlanckFilter( Settings* settings );
    virtual ~FokkerPlanckFilter();
    virtual void FilterMoments( Tensor& u ) const;
    virtual void FilterMoments( Tensor& v, const Tensor& u ) const;
    virtual void FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const;
    void SetupFilter();
};

#endif    // FOKKERPLANCKFILTER_H
