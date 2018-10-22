#include "minmod.h"
#include "math.h"

Minmod::Minmod( Settings* settings, Mesh* mesh, Closure* closure ) : Limiter( settings, mesh, closure ) {}

Minmod::~Minmod() {}

double Minmod::CalculateSlope( const double& u0, const double& u1, const double& u2 ) { return ( 1 / _dx ) * minmod( u2 - u1, u1 - u0 ); }

double Minmod::minmod( const double& a, const double& b ) {
    double y;
    if( fabs( a ) < fabs( b ) && a * b > 0 )
        y = a;
    else if( fabs( b ) < fabs( a ) && a * b > 0 )
        y = b;
    else
        y = 0;
    return y;
}
