#include "vector.cpp"

template <class T> VectorSpace::Vector<T> column( VectorSpace::Matrix<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( const VectorSpace::Matrix<T>& mat, unsigned i );

namespace VectorSpace {

template <class T> class Matrix    // column major
{
  private:
    T* _data;
    unsigned _rows;
    unsigned _columns;

  public:
    Matrix();
    Matrix( unsigned rows, unsigned columns, bool skipZeroInit = false );
    Matrix( unsigned rows, unsigned columns, T init );
    Matrix( const Matrix& other );
    ~Matrix();
    void operator=( const Matrix<T>& other );
    void resize( unsigned rows, unsigned columns );

    T& operator()( unsigned i, unsigned j );
    const T& operator()( unsigned i, unsigned j ) const;

    Matrix<T> operator+( const Matrix<T>& other ) const;
    Matrix<T> operator-( const Matrix<T>& other ) const;
    Matrix<T> operator*( const Matrix<T>& other ) const;

    Matrix<T> operator*( const T& scalar ) const;
    Matrix<T> operator/( const T& scalar ) const;
    Matrix<T> operator+( const T& scalar ) const;
    Matrix<T> operator-( const T& scalar ) const;

    Vector<T> operator*( const Vector<T>& vector ) const;

    void operator+=( const Matrix<T>& other );
    void operator-=( const Matrix<T>& other );

    double* GetPointer() { return _data; }

    unsigned rows() const;
    unsigned columns() const;
    void reset();
    bool isSymmetric() const;

    friend Vector<T> column<>( Matrix<T>& mat, unsigned i );
    friend Vector<T> column<>( const Matrix<T>& mat, unsigned i );
    friend void gesv<>( Matrix<T>& A, Vector<T>& b, int* ipiv );
    friend void posv<>( Matrix<T>& A, Vector<T>& b );
};

template <class T> Matrix<T>::Matrix() : _data( nullptr ), _rows( 0 ), _columns( 0 ) {}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, bool skipZeroInit ) : _rows( rows ), _columns( columns ) {
    _data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    if( !skipZeroInit ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            for( unsigned i = 0; i < _rows; ++i ) {
                _data[j * _rows + i] = 0.0;
            }
        }
    }
}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, T init ) : _rows( rows ), _columns( columns ) {
    _data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[j * _rows + i] = init;
        }
    }
}

template <class T> Matrix<T>::Matrix( const Matrix& other ) : _rows( other._rows ), _columns( other._columns ) {
    _data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[j * _rows + i] = other._data[j * _rows + i];
        }
    }
}

template <class T> Matrix<T>::~Matrix() {
    free( _data );    // delete[] _data;
}

template <class T> void Matrix<T>::operator=( const Matrix<T>& other ) {
    if( _data == nullptr ) {
        _rows    = other._rows;
        _columns = other._columns;
        _data    = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    }

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, other._columns );
    unsigned rows    = std::min( _rows, other._rows );

    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
            ( *this )( i, j ) = other( i, j );
        }
    }
}

template <class T> T& Matrix<T>::operator()( unsigned i, unsigned j ) { return _data[j * _rows + i]; }

template <class T> const T& Matrix<T>::operator()( unsigned i, unsigned j ) const { return _data[j * _rows + i]; }

template <class T> Matrix<T> Matrix<T>::operator+( const Matrix<T>& other ) const {
    Matrix<T> res( _rows, _columns, true );

    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) + other( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator-( const Matrix<T>& other ) const {
    Matrix<T> res( _rows, _columns, true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) - other( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const Matrix<T>& other ) const {
    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned rows = std::min( _columns, other._rows );

    Matrix<T> res( this->_rows, other.columns() );
    for( unsigned j = 0; j < other._columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned k = 0; k < rows; k++ ) {
                res( i, j ) += ( *this )( i, k ) * other( k, j );
            }
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator+( const T& scalar ) const {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) + scalar;
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator-( const T& scalar ) const {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) - scalar;
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator*( const T& scalar ) const {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) * scalar;
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator/( const T& scalar ) const {
    Matrix<T> res( _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) / scalar;
        }
    }
    return res;
}

template <class T> Vector<T> Matrix<T>::operator*( const Vector<T>& vector ) const {
    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, vector.size() );
    Vector<T> res( _rows );
    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res[i] += ( *this )( i, j ) * vector[j];
        }
    }
    return res;
}

template <class T> void Matrix<T>::operator+=( const Matrix<T>& other ) {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            ( *this )( i, j ) += other( i, j );
        }
    }
}

template <class T> void Matrix<T>::operator-=( const Matrix<T>& other ) {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            ( *this )( i, j ) -= other( i, j );
        }
    }
}

template <class T> unsigned Matrix<T>::rows() const { return _rows; }

template <class T> unsigned Matrix<T>::columns() const { return _columns; }

template <class T> void Matrix<T>::reset() {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            ( *this )( i, j ) = 0.0;
        }
    }
}

template <class T> void Matrix<T>::resize( unsigned rows, unsigned columns ) {
    auto dataOld = _data;
    _data        = static_cast<T*>( malloc( rows * columns * sizeof( T ) ) );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[j * rows + i] = dataOld[j * _rows + i];
        }
    }
    delete dataOld;
    _rows    = rows;
    _columns = columns;
}

template <class T> bool Matrix<T>::isSymmetric() const {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
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
