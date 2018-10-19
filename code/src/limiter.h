#ifndef LIMITER_H
#define LIMITER_H

#include <cpptoml.h>

#include "closure.h"
#include "mesh.h"
#include "settings.h"

class Limiter
{
  private:
    Matrix SlopeInternal( const Matrix& u0, const Matrix& u1, const Matrix& u2 );
    double SlopeBoundPres( const double& u, const double& slope );

  protected:
    Settings* _settings;
    Mesh* _mesh;
    Closure* _closure;

    double _dx;

  public:
    Limiter( Settings* settings, Mesh* mesh, Closure* closure );
    virtual ~Limiter();
    static Limiter* Create( Settings* settings, Mesh* mesh, Closure* closure );
    Matrix Slope( const Matrix& lambda1, const Matrix& lambda2, const Matrix& lambda3 );
    virtual double CalculateSlope( const double& u0, const double& u1, const double& u2 ) = 0;

  private:
    Limiter() {}
};

#endif    // LIMITER_H
