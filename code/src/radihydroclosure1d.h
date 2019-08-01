#ifndef RADIHYDROCLOSURE1D_H
#define RADIHYDROCLOSURE1D_H

#include "closure.h"

class RadiHydroClosure1D : public Closure
{
  private:
    double _gamma;
    unsigned _nMoments;
    unsigned _nHydroStates;
    RadiHydroClosure1D() = delete;

  public:
    RadiHydroClosure1D( Settings* settings );
    virtual ~RadiHydroClosure1D();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Vector& out, const Vector& Lambda, bool dummy );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void DS( Vector& ds, const Vector& u ) const;
    virtual void SolveClosure( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel );
    virtual void AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const;
    virtual void SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const;
};

#endif    // RADIHYDROCLOSURE1D_H
