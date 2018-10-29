#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

namespace VectorSpace {

template<class T>
class Matrix
{
private:
    unsigned _rows;
    unsigned _columns;
    T** _data;

 public:
    Matrix();
    Matrix(unsigned rows, unsigned columns, bool skipZeroInit = false);
    Matrix(unsigned rows, unsigned columns, T init);
    Matrix(const Matrix& other);
    ~Matrix();

    T& operator()(unsigned i, unsigned j);
    const T& operator()(unsigned i, unsigned j) const;
    Matrix operator+(const Matrix& other);
    Matrix operator-(const Matrix& other);
    Matrix operator*(const Matrix& other);
    Matrix operator*(const double& scalar);
    Vector<T> operator*(const Vector<T>& vector);
    void operator=(const Matrix<T>& other);

    unsigned rows() const;
    unsigned columns() const;
    void reset();
};

}

#endif // MATRIX_H
