#ifndef LASSOFILTER_H
#define LASSOFILTER_H

#include "filter.h"

class LassoFilter : public Filter
{
  private:
    LassoFilter() = delete;
    Vector _filterParam;
    Vector _l1Norms;

  public:
    LassoFilter( Settings* settings );
    virtual ~LassoFilter();

    virtual void SetupFilter();

    virtual void FilterMoments( Tensor& u ) const;
    virtual void FilterMoments( Tensor& v, const Tensor& u ) const;
    virtual void FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const;
};

#endif    // LASSOFILTER_H
