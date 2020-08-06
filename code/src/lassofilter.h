#ifndef LASSOFILTER_H
#define LASSOFILTER_H

#include "filter.h"

class LassoFilter : public Filter
{
  private:
    LassoFilter() = delete;
    double _lambda;
    Vector _filterParam;
    Vector _l1Norms;

  public:
    LassoFilter( Settings* settings );
    virtual ~LassoFilter();

    virtual void FilterMoments( Tensor& u ) const;
    virtual void FilterMoments( Tensor& v, const Tensor& u ) const;
    virtual double FilterMoments( double u, unsigned i ) const;
};

#endif    // LASSOFILTER_H
