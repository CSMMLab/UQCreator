#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

//#include <blaze/math/DynamicMatrix.h>
//#include <blaze/math/DynamicVector.h>

// typedef blaze::DynamicMatrix<double> Matrix;
// typedef blaze::DynamicVector<double> Vector;
// typedef blaze::DynamicVector<unsigned> VectorU;
// typedef blaze::DynamicVector<blaze::DynamicMatrix<double>> MatVec;

#include "vectorspace.cpp"

typedef VectorSpace::Matrix<double> Matrix;
typedef VectorSpace::Vector<double> Vector;
typedef VectorSpace::Vector<unsigned> VectorU;
typedef std::vector<VectorSpace::Matrix<double>> MatVec;

#endif    // TYPEDEFS_H
