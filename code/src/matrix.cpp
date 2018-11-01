#include "vector.cpp"

template <class T> VectorSpace::Vector<T> column( VectorSpace::Matrix<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( const VectorSpace::Matrix<T>& mat, unsigned i );

namespace VectorSpace {

template <class T> class Matrix    // column major
{
  private:
    unsigned _rows;
    unsigned _columns;
    T* _data;

  public:
    Matrix();
    Matrix( unsigned rows, unsigned columns, bool skipZeroInit = false );
    Matrix( unsigned rows, unsigned columns, T init );
    Matrix( const Matrix& other );
    ~Matrix();

    T& operator()( unsigned i, unsigned j );
    const T& operator()( unsigned i, unsigned j ) const;
    Matrix<T> operator+( const Matrix<T>& other ) const;
    Matrix<T> operator-( const Matrix<T>& other ) const;
    Matrix<T> operator*( const Matrix<T>& other ) const;
    Matrix<T> operator*( const double& scalar ) const;
    Vector<T> operator*( const Vector<T>& vector ) const;
    void operator=( const Matrix<T>& other );
    void operator+=( const Matrix<T>& other );
    void operator-=( const Matrix<T>& other );

    unsigned rows() const;
    unsigned columns() const;
    void reset();
    bool isSymmetric() const;

    friend Vector<T> column<>( Matrix<T>& mat, unsigned i );
    friend Vector<T> column<>( const Matrix<T>& mat, unsigned i );
    friend void gesv<>( Matrix<T>& A, Vector<T>& b, int* ipiv );
};

template <class T> Matrix<T>::Matrix() : _data( nullptr ), _rows( 0 ), _columns( 0 ) {}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, bool skipZeroInit ) : _rows( rows ), _columns( columns ) {
    _data = new T[_rows * _columns];
    if( !skipZeroInit ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned j = 0; j < _columns; ++j ) {
                _data[i * _columns + j] = 0.0;
            }
        }
    }
}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, T init ) : _rows( rows ), _columns( columns ) {
    _data = new T[_rows * _columns];
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            _data[i * _columns + j] = init;
        }
    }
}

template <class T> Matrix<T>::Matrix( const Matrix& other ) {
    _rows    = other._rows;
    _columns = other._columns;
    _data    = new T[_rows * _columns];
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            _data[i * _columns + j] = other._data[i * _columns + j];
        }
    }
}

template <class T> Matrix<T>::~Matrix() { delete[] _data; }

template <class T> T& Matrix<T>::operator()( unsigned i, unsigned j ) { return _data[i * _columns + j]; }

template <class T> const T& Matrix<T>::operator()( unsigned i, unsigned j ) const { return _data[i * _columns + j]; }

template <class T> Matrix<T> Matrix<T>::operator+( const Matrix<T>& other ) const {
    Matrix<T> res( _rows, _columns, true );
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            res( i, j ) = ( *this )( i, j ) + other( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator-( const Matrix<T>& other ) const {
    Matrix<T> res( _rows, _columns, true );
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            res( i, j ) = ( *this )( i, j ) - other( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const Matrix<T>& other ) const {
    Matrix<T> res( this->_rows, other.columns() );
    for( unsigned i = 0; i < _rows; i++ ) {
        for( unsigned j = 0; j < other.columns(); j++ ) {
            for( unsigned k = 0; k < other.rows(); k++ ) {
                res( i, j ) += ( *this )( i, k ) * other( k, j );
            }
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const double& scalar ) const {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            res( i, j ) = ( *this )( i, j ) * scalar;
        }
    }
    return res;
}

template <class T> Vector<T> Matrix<T>::operator*( const Vector<T>& vector ) const {
    Vector<T> res( _rows );
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            res[i] += ( *this )( i, j ) * vector[j];
        }
    }
    return res;
}

template <class T> void Matrix<T>::operator=( const Matrix<T>& other ) {
    if( _data == nullptr ) {
        _rows    = other._rows;
        _columns = other._columns;
        _data    = new T[_rows * _columns];
    }
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            ( *this )( i, j ) = other( i, j );
        }
    }
}

template <class T> void Matrix<T>::operator+=( const Matrix<T>& other ) {
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            ( *this )( i, j ) += other( i, j );
        }
    }
}

template <class T> void Matrix<T>::operator-=( const Matrix<T>& other ) {
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            ( *this )( i, j ) -= other( i, j );
        }
    }
}

template <class T> unsigned Matrix<T>::rows() const { return _rows; }

template <class T> unsigned Matrix<T>::columns() const { return _columns; }

template <class T> void Matrix<T>::reset() {
    for( unsigned i = 0; i < _rows; ++i ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            ( *this )( i, j ) = 0.0;
        }
    }
}

template <class T> bool Matrix<T>::isSymmetric() const {
    for( unsigned i = 1; i < _rows; ++i ) {
        for( unsigned j = i - 1; j < _columns; ++j ) {
            if( ( *this )( i, j ) != ( *this )( j, i ) ) {
                return false;
            }
        }
    }
    return true;
}

template <class T> class IdentityMatrix : public Matrix<T>
{
  private:
  public:
    IdentityMatrix();
    IdentityMatrix( unsigned n );
};

template <class T> IdentityMatrix<T>::IdentityMatrix() : Matrix<T>() {}

template <class T> IdentityMatrix<T>::IdentityMatrix( unsigned n ) : Matrix<T>( n, n, false ) {
    for( unsigned i = 0; i < n; ++i ) {
        ( *this )( i, i ) = 1;
    }
}

}    // namespace VectorSpace
