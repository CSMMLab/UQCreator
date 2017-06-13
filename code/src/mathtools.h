#ifndef MATHTOOLS_H
#define MATHTOOLS_H

#include <blaze/Blaze.h>

typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;


class MathTools
{
private:
    MathTools(){}
public:
    ~MathTools(){}
    static double pythag(const double a, const double b);
    static void computeEigenValTriDiagMatrix(const matrix CM);

};

#endif // MATHTOOLS_H
