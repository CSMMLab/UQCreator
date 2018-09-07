#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

#include <iostream>

#define PI 3.14159265359

class MathTools
{
  private:
    MathTools() {}

  public:
    ~MathTools() {}
    static double Pythag( const double a, const double b );
    static std::pair<blaze::DynamicVector<double>, blaze::DynamicMatrix<double>>
    ComputeEigenValTriDiagMatrix( const blaze::DynamicMatrix<double> CM );
};

#endif    // MATHTOOLS_H
