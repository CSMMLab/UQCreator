#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

#include "../src/vectorspace.cpp"

typedef VectorSpace::Tensor<double> Tensor;
typedef VectorSpace::Matrix<double> Matrix;
typedef VectorSpace::Vector<double> Vector;
typedef VectorSpace::Vector<unsigned> VectorU;
typedef std::vector<std::vector<unsigned>> MatrixU;
typedef std::vector<VectorSpace::Matrix<double>> MatVec;
typedef std::vector<VectorSpace::Tensor<double>> MatTens;

#endif    // TYPEDEFS_H
