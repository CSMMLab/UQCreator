#ifndef THERMALRADIATIONCLOSURE_H
#define THERMALRADIATIONCLOSURE_H

#include "closure.h"

class ThermalRadiationClosure : public Closure
{
  private:
    ThermalRadiationClosure() = delete;
    unsigned _nHydroStates;
    unsigned _nMoments;

  public:
    ThermalRadiationClosure( Settings* settings );
    virtual ~ThermalRadiationClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Vector& out, const Vector& Lambda, bool dummy );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
    virtual void SolveClosure( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel );
    virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel );
    virtual void AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const;
    virtual void SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const;
    virtual double CalcNorm( Vector& test, unsigned nTotal ) const;
};


#endif // THERMALRADIATIONCLOSURE_H
