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
    static std::pair<Vector, Matrix> ComputeEigenValTriDiagMatrix( const Matrix& CM );
    static unsigned BinomialCoefficient( unsigned n, unsigned k );
    static int Factorial( int i );
    static double csch( const double x );
    static double coth( const double x );
    static double max( double a, double b );
    static double min( double a, double b );
    unsigned max( unsigned a, unsigned b );
    static double sign( double a );
};

#endif    // MATHTOOLS_H
