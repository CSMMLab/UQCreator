#ifndef M1EQUATIONS1D_H
#define M1EQUATIONS1D_H

#include "problem.h"
#include <cmath>

class M1Equations1D : public Problem
{
  protected:
    double _sigmaA;    // absorption coefficient
    double _sigmaS;    // scattering coefficient
    double _sigmaT;    // total crossection

    double ComputeAlpha( const double u1Du0 ) const;
    double RootFun( const double alpha, const double u1Du0 ) const;
    double Bisection( double alphaA, double alphaB, const double u1Du0 ) const;

  public:
    /**
     * @brief M1Equations1D constructur without setting up system matrices used for radiation hydrodynamics
     * @param settings Settings pointer
     * @param noSystemMatrix dummy bool
     */
    M1Equations1D( Settings* settings );
    virtual ~M1Equations1D();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif    // M1EQUATIONS1D_H
