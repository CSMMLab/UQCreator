#ifndef EULER2DFPFILTER_H
#define EULER2DFPFILTER_H

#include "eulerclosure2d.h"

class Euler2DFPFilter : public EulerClosure2D
{
  private:
    double _lambda;
    Vector _filterFunction;
    Euler2DFPFilter() = delete;
    Matrix Filter( const Matrix& u, unsigned nTotal ) const;

  public:
    Euler2DFPFilter( Settings* settings );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
    // virtual ~Euler2DFPFilter();
};

#endif    // EULER2DFPFILTER_H
