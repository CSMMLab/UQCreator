#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

#include "vectorspace.cpp"

typedef VectorSpace::Matrix<double> Matrix;
typedef VectorSpace::Vector<double> Vector;
typedef VectorSpace::Vector<unsigned> VectorU;
typedef std::vector<std::vector<unsigned>> MatrixU;
typedef std::vector<VectorSpace::Matrix<double>> MatVec;

#endif    // TYPEDEFS_H
