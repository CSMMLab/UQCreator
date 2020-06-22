#include "vector.cpp"

template <class T> VectorSpace::Vector<T> column( VectorSpace::Tensor<T>& mat, unsigned i );
template <class T> VectorSpace::Vector<T> column( const VectorSpace::Tensor<T>& mat, unsigned i );

namespace VectorSpace {

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
    void resize( unsigned frontRows, unsigned rows, unsigned columns );
    Tensor<T> Add( const Tensor<T>& other, unsigned frontRows, unsigned rows, unsigned columns ) const;

    T& operator()( unsigned l, unsigned i, unsigned j );
    const T& operator()( unsigned l, unsigned i, unsigned j ) const;

    Tensor<T> operator+( const Tensor<T>& other ) const;
    Tensor<T> operator-( const Tensor<T>& other ) const;
    // Tensor<T> operator*( const Tensor<T>& other ) const;

    Tensor<T> operator*( const T& scalar ) const;
    Tensor<T> operator/( const T& scalar ) const;
    Tensor<T> operator+( const T& scalar ) const;
    Tensor<T> operator-( const T& scalar ) const;

    // Vector<T> operator*( const Vector<T>& vector ) const;

    void operator+=( const Tensor<T>& other );
    void operator-=( const Tensor<T>& other );

    double* GetPointer() { return _data; }

    unsigned rows() const;
    unsigned columns() const;
    void reset();

    friend Vector<T> column<>( Tensor<T>& mat, unsigned i );
    friend Vector<T> column<>( const Tensor<T>& mat, unsigned i );
};

template <class T> Tensor<T>::Tensor() : _data( nullptr ), _rows( 0 ), _columns( 0 ), _frontRows( 0 ) {}

template <class T>
Tensor<T>::Tensor( unsigned frontRows, unsigned rows, unsigned columns, bool skipZeroInit )
    : _rows( rows ), _columns( columns ), _frontRows( frontRows ) {
    _data = static_cast<T*>( malloc( _frontRows * _rows * _columns * sizeof( T ) ) );
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
Tensor<T>::Tensor( unsigned frontRows, unsigned rows, unsigned columns, T init ) : _rows( rows ), _columns( columns ), _frontRows( frontRows ) {
    _data = static_cast<T*>( malloc( _frontRows * _rows * _columns * sizeof( T ) ) );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                _data[( j * _rows + i ) * _frontRows + l] = init;
            }
        }
    }
}

template <class T> Tensor<T>::Tensor( const Tensor& other ) : _rows( other._rows ), _columns( other._columns ), _frontRows( other._frontRows ) {
    _data = static_cast<T*>( malloc( _frontRows * _rows * _columns * sizeof( T ) ) );
    for( unsigned j = 0; j < _columns; ++j ) {
        for( unsigned i = 0; i < _rows; ++i ) {
            for( unsigned l = 0; l < _frontRows; ++l ) {
                _data[( j * _rows + i ) * _frontRows + l] = other._data[( j * _rows + i ) * _frontRows + l];
            }
        }
    }
}

template <class T> Tensor<T>::~Tensor() {
    free( _data );    // delete[] _data;
}

template <class T> void Tensor<T>::operator=( const Tensor<T>& other ) {
    if( _data == nullptr ) {
        _rows      = other._rows;
        _columns   = other._columns;
        _frontRows = other._frontRows;
        _data      = static_cast<T*>( malloc( _rows * _columns * sizeof( T ) ) );
    }

    // ensure that only allocated data is written when sizes differ (due to adaptivity)
    unsigned columns    = std::min( _columns, other._columns );
    unsigned rows       = std::min( _rows, other._rows );
    unsigned _frontRows = std::min( _frontRows, other._frontRows );

    for( unsigned l = 0; l < _frontRows; ++l ) {
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

template <class T> T& Tensor<T>::operator()( unsigned l, unsigned i, unsigned j ) { _data[( j * _rows + i ) * _frontRows + l]; }

template <class T> const T& Tensor<T>::operator()( unsigned l, unsigned i, unsigned j ) const { return _data[( j * _rows + i ) * _frontRows + l]; }

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
/*
template <class T> void Tensor<T>::resize( unsigned rows, unsigned columns ) {
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
}*/

}    // namespace VectorSpace
