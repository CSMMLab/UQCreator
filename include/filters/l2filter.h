#ifndef L2FILTER_H
#define L2FILTER_H

#include "filter.h"

class L2Filter : public Filter
{
  private:
    L2Filter() = delete;

  public:
    L2Filter( Settings* settings );
    virtual ~L2Filter();
    void SetupFilter();

    virtual void FilterMoments( Tensor& u ) const;
    virtual void FilterMoments( Tensor& v, const Tensor& u ) const;
    virtual void FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const;
};

#endif    // L2FILTER_H
