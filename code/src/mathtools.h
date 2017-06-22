#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <blaze/Blaze.h>

#define PI 3.14159265359

typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;


class MathTools
{
private:
    MathTools(){}
public:
    ~MathTools(){}
    static double pythag(const double a, const double b);
    static std::pair<vector,matrix> computeEigenValTriDiagMatrix(const matrix CM);

};

#endif // MATHTOOLS_H
