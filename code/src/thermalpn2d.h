#ifndef THERMALPN2D_H
#define THERMALPN2D_H

#include "mathtools.h"
#include "problem.h"

class ThermalPN2D : public Problem
{
  private:
    double _epsilon;
    double _c;
    double _cV;
    double _a;
    double _TRef;
    double _sigma;
    double _alpha;
    bool _suOlson;
    unsigned _constitutiveLaw;
    Matrix _AbsAx;
    Matrix _AbsAz;
    std::vector<Vector> _xiQuad;
    std::vector<double> _variances;
    double ScaledInternalEnergy( double TTilde ) const;
    double ScaledTemperature( double eTilde ) const;
    double Delta( int l, int k ) const;
    void SetupSystemMatrices();
    int GlobalIndex( int l, int k ) const;
    Matrix _Ax;    // flux matrices
    Matrix _Ay;
    Matrix _Az;
    unsigned _nMoments;    // number of radiative moments
    unsigned _N;

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

    int Sgn( int k ) const;
    int kPlus( int k ) const;
    int kMinus( int k ) const;

  public:
    ThermalPN2D( Settings* settings );
    virtual ~ThermalPN2D();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif // THERMALPN2D_H
