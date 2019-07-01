#ifndef RADIATIONHYDRODYNAMICS_H
#define RADIATIONHYDRODYNAMICS_H

#include "pnequations.h"
#include "mathtools.h"
#include <cmath>

class RadiationHydrodynamics : public PNEquations
{
    double _c; // radiation velocity
    double _P;
    double _R;
    double _gamma;
    unsigned _nMoments;
    double Delta( int l, int k ) const;

    void SetupSystemMatrices();
    virtual Matrix Source( const Matrix& uQ ) const;
    double Er0(const Vector& u)const;
    Vector Fr0(const Vector& u)const;
    double SE(const Vector& u)const;
    Vector SF(const Vector& u)const;
public:
    RadiationHydrodynamics(Settings* settings );
    Matrix F( const Vector& u );
};

#endif // RADIATIONHYDRODYNAMICS_H
