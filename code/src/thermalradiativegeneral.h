#ifndef THERMALRADIATIVEGENERAL_H
#define THERMALRADIATIVEGENERAL_H

#include "problem.h"

class ThermalRadiativeGeneral : public Problem
{
  private:
    double _epsilon;
    double _c;
    double _a;
    double _TRef;
    double _sigma;
    double _alpha;
    bool _suOlson;
    std::vector<Vector> _xiQuad;
    std::vector<double> _variances;
    double ScaledInternalEnergy(double TTilde )const;
    double ScaledTemperature(double eTilde ) const;

  public:
    ThermalRadiativeGeneral( Settings* settings );
    virtual ~ThermalRadiativeGeneral();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
    virtual Vector LoadIC( const Vector& x, const Vector& xi );
};

#endif // THERMALRADIATIVEGENERAL_H
