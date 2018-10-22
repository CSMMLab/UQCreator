#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <iostream>

#include "typedefs.h"

#define PI 3.14159265359

class MathTools
{
  private:
    MathTools() {}

  public:
    ~MathTools() {}
    static double Pythag( const double a, const double b );
    static std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix CM );
};

#endif    // MATHTOOLS_H
