#include "matrix.h"

namespace VectorSpace {

template <class T> Matrix<T>::Matrix() : _rows( 0 ), _columns( 0 ) {}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, bool skipZeroInit ) : _rows( rows ), _columns( columns ) {
    _data = new T[_rows];
    for( unsigned i = 0; i < _columns; ++i ) {
        _data[i] = new T[_columns];
    }
    if( !skipZeroInit ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            for( unsigned i = 0; i < _rows; ++i ) {
                _data[i][j] = 0.0;
            }
        }
    }
}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, T init ) : _rows( rows ), _columns( columns ) {
    _data = new T[_rows];
    for( unsigned i = 0; i < _columns; ++i ) {
        _data[i] = new T[_columns];
    }
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[i][j] = init;
        }
    }
}

template <class T> Matrix<T>::~Matrix() {
    for( unsigned i = 0; i < _rows; ++i ) {
        delete[] _data[i];
    }
    delete[] _data;
}

template <class T> Matrix<T>::Matrix( const Matrix& other ) {
    _rows    = other._rows;
    _columns = other._columns;
    _data    = new T[_rows];
    for( unsigned i = 0; i < _columns; ++i ) {
        _data[i] = new T[_columns];
    }
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[i][j] = other._data[i][j];
        }
    }
}

template <class T> T& Matrix<T>::operator()( unsigned i, unsigned j ) { return _data[i][j]; }

template <class T> const T& Matrix<T>::operator()( unsigned i, unsigned j ) const { return _data[i][j]; }

template <class T> Matrix<T> Matrix<T>::operator+( const Matrix<T>& other ) {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res = this->_data[i][j] + other._data[i][j];
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator-( const Matrix<T>& other ) {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res = this->_data[i][j] - other._data[i][j];
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const Matrix<T>& other ) {
    Matrix<T> res( this->rows, other.columns() );
    for( unsigned i = 0; i < _rows; i++ ) {
        for( unsigned j = 0; j < other.columns(); j++ ) {
            for( unsigned k = 0; k < other.rows(); k++ ) {
                res( i, j ) += this( i, k ) * other( k, j );
            }
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const double& scalar ) {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res = this->_data[i][j] * scalar;
        }
    }
    return res;
}

template <class T> Vector<T> Matrix<T>::operator*( const Vector<T>& vector ) {
    Vector<T> res( _rows, true );
    for( unsigned i = 0; i < _rows; ++i ) {
        res[i] = dot( column( this, i ), vector );
    }
    return res;
}

template <class T> void Matrix<T>::operator=( const Matrix<T>& other ) {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            this->_data[i][j] = other._data[i][j];
        }
    }
}

template <class T> unsigned Matrix<T>::rows() const { return _rows; }

template <class T> unsigned Matrix<T>::columns() const { return _columns; }

template <class T> void Matrix<T>::reset() {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[i][j] = 0.0;
        }
    }
}

}    // namespace VectorSpace
