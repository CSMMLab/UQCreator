#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

#include "vectorspace.cpp"

typedef VectorSpace::Matrix<double> Matrix;
typedef VectorSpace::Vector<double> Vector;
typedef VectorSpace::Vector<unsigned> VectorU;
typedef std::vector<VectorSpace::Matrix<double>> MatVec;
typedef VectorSpace::FluxMatrix<double> FluxMatrix;

#endif    // TYPEDEFS_H
