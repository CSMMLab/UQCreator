#ifndef RADIATIONHYDRODYNAMICS_H
#define RADIATIONHYDRODYNAMICS_H

#include "mathtools.h"
#include "pnequations.h"
#include <cmath>

class RadiationHydrodynamics : public PNEquations
{
    double _c;    // radiation velocity
    double _P;
    double _R;
    double _gamma;
    unsigned _nMoments;
    double Delta( int l, int k ) const;

    void SetupSystemMatrices();
    double Er0( const Vector& u ) const;
    Vector Fr0( const Vector& u ) const;
    double SE( const Vector& u ) const;
    Vector SF( const Vector& u ) const;
    Matrix FRadiation( const Vector& u );
    Matrix FEuler( const Vector& u );
    Matrix F( const Vector& u );

  public:
    RadiationHydrodynamics( Settings* settings );

    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level );
    virtual Matrix Source( const Matrix& uQ ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
};

#endif    // RADIATIONHYDRODYNAMICS_H
