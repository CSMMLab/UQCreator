#ifndef EULER1DFPFILTER_H
#define EULER1DFPFILTER_H

#include "eulerclosure.h"

class Euler1DFPFilter : public EulerClosure
{
  private:
    double _lambda;
    Vector _filterFunction;
    Euler1DFPFilter() = delete;
    Matrix Filter( const Matrix& u, unsigned nTotal ) const;
    bool _hyperbolicityFix;

  public:
    Euler1DFPFilter( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
    // virtual ~Euler1DFPFilter();
};

#endif // EULER1DFPFILTER_H
