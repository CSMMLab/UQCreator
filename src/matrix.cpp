#include "vector.cpp"
#include <assert.h>

template <class T> VectorSpace::Vector<T> column( VectorSpace::Matrix<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( const VectorSpace::Matrix<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( VectorSpace::Tensor<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( const VectorSpace::Tensor<T>& mat, unsigned i );

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
    Matrix<T> Add( const Matrix<T>& other, unsigned rows, unsigned columns ) const;

    T& operator()( unsigned i, unsigned j );
    const T& operator()( unsigned i, unsigned j ) const;

    Matrix<T> operator+( const Matrix<T>& other ) const;
    Matrix<T> operator-( const Matrix<T>& other ) const;
    Matrix<T> operator*( const Matrix<T>& other ) const;

    Matrix<T> operator*( const T& scalar ) const;
    Matrix<T> operator/( const T& scalar ) const;
    Matrix<T> operator+( const T& scalar ) const;
    Matrix<T> operator-( const T& scalar ) const;
    Matrix<T> transpose() const;
    Matrix<T> inv() const;

    Vector<T> operator*( const Vector<T>& vector ) const;

    void operator+=( const Matrix<T>& other );
    void operator-=( const Matrix<T>& other );

    double* GetPointer() { return _data; }
    double* GetPointer() const { return _data; }

    unsigned rows() const;
    unsigned columns() const;
    void reset();

    bool isSymmetric() const;

    friend Vector<T> column<>( Matrix<T>& mat, unsigned i );
    friend Vector<T> column<>( const Matrix<T>& mat, unsigned i );
    friend void gesv<>( const Matrix<T>& A, Vector<T>& b, int* ipiv );
    friend void posv<>( Matrix<T>& A, Vector<T>& b );
};

template <class T> Matrix<T>::Matrix() : _data( nullptr ), _rows( 0 ), _columns( 0 ) {}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, bool skipZeroInit ) : _data( nullptr ), _rows( rows ), _columns( columns ) {
    //_data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    _data = new T[_rows * _columns];
    if( !skipZeroInit ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            for( unsigned i = 0; i < _rows; ++i ) {
                assert( j * _rows + i < _rows * _columns );
                _data[j * _rows + i] = 0.0;
            }
        }
    }
}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, T init ) : _data( nullptr ), _rows( rows ), _columns( columns ) {
    //_data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    _data = new T[_rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            assert( j * _rows + i < _rows * _columns );
            _data[j * _rows + i] = init;
        }
    }
}

template <class T> Matrix<T>::Matrix( const Matrix& other ) : _data( nullptr ), _rows( other._rows ), _columns( other._columns ) {
    //_data = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    _data = new T[_rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            assert( j * _rows + i < _rows * _columns );
            _data[j * _rows + i] = other._data[j * _rows + i];
        }
    }
}

template <class T> Matrix<T>::~Matrix() { delete[] _data; }

template <class T> void Matrix<T>::operator=( const Matrix<T>& other ) {
    if( _data == nullptr ) {
        _rows    = other._rows;
        _columns = other._columns;
        _data    = new T[_rows * _columns];
    }

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, other._columns );
    unsigned rows    = std::min( _rows, other._rows );

    for( unsigned i = 0; i < rows; ++i ) {
        for( unsigned j = 0; j < columns; ++j ) {
            ( *this )( i, j ) = other( i, j );
        }
    }
}

template <class T> Matrix<T> Matrix<T>::Add( const Matrix<T>& other, unsigned rows, unsigned columns ) const {
    Matrix<T> res( _rows, _columns, true );
    assert( _rows == other._rows );
    assert( _columns == other.columns() );

    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) + other( i, j );
        }
    }
    return res;
}

template <class T> T& Matrix<T>::operator()( unsigned i, unsigned j ) {
#ifndef NDEBUG
    assert( i < _rows );
    assert( j < _columns );
    assert( j * _rows + i < _rows * _columns );
#endif
    return _data[j * _rows + i];
}

template <class T> const T& Matrix<T>::operator()( unsigned i, unsigned j ) const {
#ifndef NDEBUG
    assert( i < _rows );
    assert( j < _columns );
    assert( j * _rows + i < _rows * _columns );
#endif
    return _data[j * _rows + i];
}

template <class T> Matrix<T> Matrix<T>::operator+( const Matrix<T>& other ) const {

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, other._columns );
    unsigned rows    = std::min( _rows, other._rows );

    Matrix<T> res( rows, columns, true );

    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
            res( i, j ) = ( *this )( i, j ) + other( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::operator-( const Matrix<T>& other ) const {

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, other._columns );
    unsigned rows    = std::min( _rows, other._rows );

    Matrix<T> res( rows, columns, true );

    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
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

template <class T> Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> res( _columns, _rows, true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            res( j, i ) = ( *this )( i, j );
        }
    }
    return res;
}

template <class T> Matrix<T> Matrix<T>::inv() const {
    Matrix<T> res( _rows, _columns, true );
    Vector<T> v( _rows );
    int* ipiv = new int[_rows];
    for( unsigned i = 0; i < _columns; ++i ) {
        v.reset();
        v[i] = 1.0;
        gesv( *this, v, ipiv );
        for( unsigned j = 0; j < _columns; ++j ) res( j, i ) = v[j];
    }
    delete[] ipiv;
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
    unsigned columns = std::min( _columns, other._columns );
    unsigned rows    = std::min( _rows, other._rows );
    // assert( _rows == other._rows );
    // assert( _columns == other.columns() );
    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
            ( *this )( i, j ) += other( i, j );
        }
    }
}

template <class T> void Matrix<T>::operator-=( const Matrix<T>& other ) {
    assert( _rows == other._rows );
    assert( _columns == other.columns() );
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
    T* dataOld = _data;
    //_data        = static_cast<T*>( malloc( rows * columns * sizeof( T ) ) );
    _data = new T[_rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[j * rows + i] = dataOld[j * _rows + i];
        }
    }
    delete[] dataOld;
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

template <class T> class Tensor    // column major
{
  private:
    T* _data;
    unsigned _rows;
    unsigned _columns;
    unsigned _frontRows;

  public:
    Tensor();
    Tensor( unsigned frontRows, unsigned rows, unsigned columns, bool skipZeroInit = false );
    Tensor( unsigned frontRows, unsigned rows, unsigned columns, T init );
    Tensor( const Tensor& other );
    ~Tensor();
    void operator=( const Tensor<T>& other );
    void resize( unsigned frontrows, unsigned rows, unsigned columns );
    Tensor<T> Add( const Tensor<T>& other, unsigned frontRows, unsigned rows, unsigned columns ) const;

    T& operator()( unsigned l, unsigned i, unsigned j );
    const T& operator()( unsigned l, unsigned i, unsigned j ) const;

    Tensor<T> operator+( const Tensor<T>& other ) const;
    Tensor<T> operator-( const Tensor<T>& other ) const;
    // Tensor<T> operator*( const Tensor<T>& other ) const;
    Tensor<T> operator*( const Matrix<T>& other ) const;

    Tensor<T> operator*( const T& scalar ) const;
    Tensor<T> operator/( const T& scalar ) const;
    Tensor<T> operator+( const T& scalar ) const;
    Tensor<T> operator-( const T& scalar ) const;

    // Vector<T> operator*( const Vector<T>& vector ) const;

    void operator+=( const Tensor<T>& other );
    void operator-=( const Tensor<T>& other );

    double* GetPointer() { return _data; }
    double* GetPointer() const { return _data; }

    unsigned frontRows() const;
    unsigned rows() const;
    unsigned columns() const;
    void reset();
    // reset tensor starting, leaving out indices 0,...,startFRow-1; ...
    void reset( unsigned startFRow, unsigned startRow, unsigned startColumns );

    friend Vector<T> column<>( Tensor<T>& mat, unsigned i );
    friend Vector<T> column<>( const Tensor<T>& mat, unsigned i );
};

template <class T> Tensor<T>::Tensor() : _data( nullptr ), _rows( 0 ), _columns( 0 ), _frontRows( 0 ) {}

template <class T>
Tensor<T>::Tensor( unsigned frontRows, unsigned rows, unsigned columns, bool skipZeroInit )
    : _data( nullptr ), _rows( rows ), _columns( columns ), _frontRows( frontRows ) {
    _data = new T[_frontRows * _rows * _columns];
    if( !skipZeroInit ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            for( unsigned i = 0; i < _rows; ++i ) {
                for( unsigned l = 0; l < _frontRows; ++l ) {
                    _data[( j * _rows + i ) * _frontRows + l] = 0.0;
                }
            }
        }
    }
}

template <class T>
Tensor<T>::Tensor( unsigned frontRows, unsigned rows, unsigned columns, T init )
    : _data( nullptr ), _rows( rows ), _columns( columns ), _frontRows( frontRows ) {
    _data = new T[_frontRows * _rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                _data[( j * _rows + i ) * _frontRows + l] = init;
            }
        }
    }
}

template <class T>
Tensor<T>::Tensor( const Tensor& other ) : _data( nullptr ), _rows( other._rows ), _columns( other._columns ), _frontRows( other._frontRows ) {
    _data = new T[_frontRows * _rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                _data[( j * _rows + i ) * _frontRows + l] = other._data[( j * _rows + i ) * _frontRows + l];
            }
        }
    }
}

template <class T> Tensor<T>::~Tensor() {
    if( _data ) delete[] _data;
}

template <class T> void Tensor<T>::operator=( const Tensor<T>& other ) {
    if( _data == nullptr ) {
        _rows      = other._rows;
        _columns   = other._columns;
        _frontRows = other._frontRows;
        //_data      = static_cast<T*>( malloc( _frontRows * _rows * _columns * sizeof( T ) ) );
        _data = new T[_frontRows * _rows * _columns];
    }

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns   = std::min( _columns, other._columns );
    unsigned rows      = std::min( _rows, other._rows );
    unsigned frontRows = std::min( _frontRows, other._frontRows );

    for( unsigned l = 0; l < frontRows; ++l ) {
        for( unsigned i = 0; i < rows; ++i ) {
            for( unsigned j = 0; j < columns; ++j ) {
                ( *this )( l, i, j ) = other( l, i, j );
            }
        }
    }
}

template <class T> Tensor<T> Tensor<T>::Add( const Tensor<T>& other, unsigned frontRows, unsigned rows, unsigned columns ) const {
    Tensor<T> res( _frontRows, _rows, _columns, true );

    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < rows; ++i ) {
            for( unsigned l = 0; l < frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) + other( l, i, j );
            }
        }
    }
    return res;
}

template <class T> T& Tensor<T>::operator()( unsigned l, unsigned i, unsigned j ) {
    assert( ( j * _rows + i ) * _frontRows + l < _frontRows * _rows * _columns );
    return _data[( j * _rows + i ) * _frontRows + l];
}

template <class T> const T& Tensor<T>::operator()( unsigned l, unsigned i, unsigned j ) const {
    assert( ( j * _rows + i ) * _frontRows + l < _frontRows * _rows * _columns );
    return _data[( j * _rows + i ) * _frontRows + l];
}

// template <class T> const T& Matrix<T>::operator()( unsigned i, unsigned j ) const { return _data[j * _rows + i]; }

template <class T> Tensor<T> Tensor<T>::operator+( const Tensor<T>& other ) const {
    Tensor<T> res( _rows, _columns, true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) + other( l, i, j );
            }
        }
    }
    return res;
}

template <class T> Tensor<T> Tensor<T>::operator-( const Tensor<T>& other ) const {
    Tensor<T> res( _rows, _columns, true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) - other( l, i, j );
            }
        }
    }
    return res;
}
/*
template <class T> Tensor<T> Tensor<T>::operator*( const Tensor<T>& other ) const {
    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned rows      = std::min( _columns, other._rows );
    unsigned frontRows = std::min( _columns, other._frontRows );

    Tensor<T> res( this->_frontRows, this->_rows, other.columns() );

    for( unsigned j = 0; j < other._columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned k = 0; k < rows; k++ ) {
                for( unsigned n = 0; n < _frontRows; ++n ) {
                    for( unsigned m = 0; m < frontRows; ++m ) {
                        res( i, j ) += ( *this )( i, k ) * other( n, k, j );
                    }
                }
            }
        }
    }
    return res;
}*/

template <class T> Tensor<T> Tensor<T>::operator+( const T& scalar ) const {
    Tensor<T> res( _frontRows, _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) + scalar;
            }
        }
    }
    return res;
}

template <class T> Tensor<T> Tensor<T>::operator-( const T& scalar ) const {
    Tensor<T> res( _frontRows, _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) - scalar;
            }
        }
    }
    return res;
}

template <class T> Tensor<T> Tensor<T>::operator*( const T& scalar ) const {
    Tensor<T> res( _frontRows, _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) * scalar;
            }
        }
    }
    return res;
}

template <class T> Tensor<T> Tensor<T>::operator*( const Matrix<T>& other ) const {
    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned rows = std::min( _columns, other.rows() );
    Tensor<T> res( _frontRows, _rows, other.columns(), 0.0 );
    for( unsigned j = 0; j < rows; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                for( unsigned n = 0; n < other.columns(); ++n ) {
                    res( l, i, n ) = res( l, i, n ) + ( *this )( l, i, j ) * other( j, n );
                }
            }
        }
    }
    return res;
}

template <class T> Tensor<T> Tensor<T>::operator/( const T& scalar ) const {
    Tensor<T> res( _frontRows, _rows, columns(), true );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res( l, i, j ) = ( *this )( l, i, j ) / scalar;
            }
        }
    }
    return res;
}
/*
template <class T> Vector<T> Tensor<T>::operator*( const Vector<T>& vector ) const {
    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns = std::min( _columns, vector.size() );
    Vector<T> res( _rows );
    for( unsigned j = 0; j < columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                res[i] += ( *this )( i, j ) * vector[j];
            }
        }
    }
    return res;
}*/

template <class T> void Tensor<T>::operator+=( const Tensor<T>& other ) {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                ( *this )( l, i, j ) += other( l, i, j );
            }
        }
    }
}

template <class T> void Tensor<T>::operator-=( const Tensor<T>& other ) {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                ( *this )( l, i, j ) -= other( l, i, j );
            }
        }
    }
}

template <class T> unsigned Tensor<T>::frontRows() const { return _frontRows; }

template <class T> unsigned Tensor<T>::rows() const { return _rows; }

template <class T> unsigned Tensor<T>::columns() const { return _columns; }

template <class T> void Tensor<T>::reset() {
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                ( *this )( l, i, j ) = 0.0;
            }
        }
    }
}

template <class T> void Tensor<T>::reset( unsigned startFRow, unsigned startRow, unsigned startColumns ) {
    for( unsigned j = startColumns; j < _columns; ++j ) {
        for( unsigned i = startRow; i < _rows; ++i ) {
            for( unsigned l = startFRow; l < _frontRows; ++l ) {
                ( *this )( l, i, j ) = 0.0;
            }
        }
    }
}

template <class T> void Tensor<T>::resize( unsigned frontrows, unsigned rows, unsigned columns ) {
    T* dataOld = _data;
    _data      = new T[rows * columns * frontrows];

    unsigned rowsMin  = std::min( rows, _rows );
    unsigned colMin   = std::min( columns, _columns );
    unsigned frontMin = std::min( frontrows, _frontRows );
    for( unsigned j = 0; j < colMin; ++j ) {
        for( unsigned i = 0; i < rowsMin; ++i ) {
            for( unsigned l = 0; l < frontMin; ++l ) {
                _data[( j * _rows + i ) * _frontRows + l] = dataOld[( j * _rows + i ) * _frontRows + l];
            }
        }
    }
    delete[] dataOld;
    _rows      = rows;
    _columns   = columns;
    _frontRows = frontrows;
}

}    // namespace VectorSpace