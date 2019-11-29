#ifndef PNEQUATIONS1D_H
#define PNEQUATIONS1D_H

#include "problem.h"
#include <cmath>

class PNEquations1D : public Problem
{
  protected:
    // moment orders for P_N
    int _N;
    double _sigmaA;    // absorption coefficient
    double _sigmaS;    // scattering coefficient
    double _sigmaT;    // total crossection

    // System Matrix for x, y and z flux
    Matrix _Ax;
    Matrix _Ay;
    Matrix _Az;

    // parameter functions for setting up system matrix
    double AParam( int l, int k ) const;
    double BParam( int l, int k ) const;
    double CParam( int l, int k ) const;
    double DParam( int l, int k ) const;
    double EParam( int l, int k ) const;
    double FParam( int l, int k ) const;

    double CTilde( int l, int k ) const;
    double DTilde( int l, int k ) const;
    double ETilde( int l, int k ) const;
    double FTilde( int l, int k ) const;

    // mathematical + index functions
    int Sgn( int k ) const;
    int kPlus( int k ) const;
    int kMinus( int k ) const;
    virtual int GlobalIndex( int l, int k ) const;

    // function for setting up system matrices
    void SetupSystemMatrices();

  public:
    PNEquations1D( Settings* settings );
    /**
     * @brief PNEquations1D constructur without setting up system matrices used for radiation hydrodynamics
     * @param settings Settings pointer
     * @param noSystemMatrix dummy bool
     */
    PNEquations1D( Settings* settings, bool noSystemMatrix );
    virtual ~PNEquations1D();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual Matrix Source( const Matrix& uQ ) const;
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif // PNEQUATIONS1D_H