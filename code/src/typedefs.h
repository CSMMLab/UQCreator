#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>

typedef blaze::DynamicMatrix<double> Matrix;
typedef blaze::DynamicVector<double> Vector;
// typedef std::vector<blaze::DynamicMatrix<double>> MatVec;
typedef blaze::DynamicVector<blaze::DynamicMatrix<double>> MatVec;

#endif    // TYPEDEFS_H
