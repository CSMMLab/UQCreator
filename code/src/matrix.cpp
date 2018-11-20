#include "vector.cpp"

#include <vector>

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
    _data = new T[_rows * _columns];
    if( !skipZeroInit ) {
        for( unsigned j = 0; j < _columns; ++j ) {
            for( unsigned i = 0; i < _rows; ++i ) {
                _data[j * _rows + i] = 0.0;
            }
        }
    }
}

template <class T> Matrix<T>::Matrix( unsigned rows, unsigned columns, T init ) : _rows( rows ), _columns( columns ) {
    _data = new T[_rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            _data[j * _rows + i] = init;
        }
    }
}

template <class T> Matrix<T>::Matrix( const Matrix& other ) : _rows( other._rows ), _columns( other._columns ) {
    _data = new T[_rows * _columns];
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
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

    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
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
    Matrix<T> res( this->_rows, other.columns() );
    for( unsigned j = 0; j < other._columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned k = 0; k < other.rows(); k++ ) {
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
    Vector<T> res( _rows );
    for( unsigned j = 0; j < _columns; ++j ) {
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

/////////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////////

template <class T> class FluxMatrix    // sparse, symmetric (upper = -lower), zero diagonal
{
  private:
    unsigned _dim;
    std::vector<T>* _data;
    std::vector<unsigned>* _columnID;
    std::vector<unsigned>* _rowID;

    void insert( T val, unsigned index, unsigned i, unsigned j );
    void remove( unsigned index, unsigned j );

    FluxMatrix() = delete;

  public:
    FluxMatrix( unsigned n );
    FluxMatrix( const FluxMatrix& other );
    ~FluxMatrix();
    void operator=( const FluxMatrix<T>& other );

    T operator()( unsigned i, unsigned j ) const;
    void set( T value, unsigned i, unsigned j );

    unsigned rows() const;
    unsigned columns() const;
    unsigned nnz() const;
};

template <class T> FluxMatrix<T>::FluxMatrix( unsigned dim ) : _dim( dim ), _data( nullptr ), _columnID( nullptr ) {
    this->_rowID = new std::vector<unsigned>( _dim + 1, 1 );
}

template <class T> FluxMatrix<T>::FluxMatrix( const FluxMatrix& other ) {
    _dim      = other._dim;
    _data     = other._data;
    _columnID = other.columnID;
    _rowID    = other.rowID;
}

template <class T> FluxMatrix<T>::~FluxMatrix() {
    if( _data != nullptr ) {
        delete _data;
        delete _columnID;
        delete _rowID;
    }
}

template <class T> void FluxMatrix<T>::set( T value, unsigned i, unsigned j ) {
    unsigned pos     = ( *_rowID )[i] - 1;
    unsigned currCol = 0;
    for( ; pos < ( *_rowID )[i] - 1; ++pos ) {
        currCol = ( *_columnID )[pos];
        if( currCol >= j ) {
            break;
        }
    }
    if( currCol != j ) {
        if( value != T() ) {
            if( i > j ) {
                std::cout << value << "\t" << pos << "\t" << i << "\t" << j << std::endl;
                this->insert( value, pos, i, j );
            }
            else if( j > i ) {
                std::cout << value << "\t" << pos << "\t" << i << "\t" << j << std::endl;
                this->insert( -value, pos, j, i );
            }
        }
    }
    else if( value == T() ) {
        this->remove( pos, i );
    }
    else {
        ( *_data )[pos] = value;
    }
}

template <class T> T FluxMatrix<T>::operator()( unsigned i, unsigned j ) const {
    std::cout << ( *_rowID )[i] - 1 << std::endl;
    for( unsigned pos = ( *_rowID )[i] - 1; pos < ( *_rowID )[i + 1] - 1; ++pos ) {
        if( ( *_columnID )[pos] == j ) {
            if( i > j )
                return ( *_data )[pos];
            else if( j > i )
                return -( *_data )[pos];
            else
                return T();
        }
        else if( ( *_columnID )[pos] > j ) {
            break;
        }
    }
    return T();    // type default
}

template <typename T> void FluxMatrix<T>::insert( T value, unsigned index, unsigned i, unsigned j ) {
    if( _data == nullptr ) {
        _data     = new std::vector<T>( 1, value );
        _columnID = new std::vector<unsigned>( 1, j );
    }
    else {
        _data->insert( _data->begin() + index, value );
        _columnID->insert( _columnID->begin() + index, j );
    }

    for( unsigned k = i; k < _dim + 1; ++k ) {
        ( *_rowID )[k] += 1;
    }
}

template <typename T> void FluxMatrix<T>::remove( unsigned index, unsigned j ) {
    _data->erase( _data->begin() + index );
    _columnID->erase( _columnID->begin() + index );

    for( unsigned k = j; k < _dim + 1; ++k ) {
        ( *_rowID )[k] -= 1;
    }
}

template <class T> unsigned FluxMatrix<T>::rows() const { return _dim; }

template <class T> unsigned FluxMatrix<T>::columns() const { return _dim; }

template <class T> unsigned FluxMatrix<T>::nnz() const { return _data->size() * 2; }

}    // namespace VectorSpace
