#include <assert.h>
#include <cmath>
#include <initializer_list>
#include <iostream>
#include <sstream>
#include <vector>

namespace VectorSpace {
template <class T> class Vector;
template <class T> class Matrix;
template <class T> class Tensor;
}    // namespace VectorSpace

template <class T> inline void gesv( const VectorSpace::Matrix<T>& A, VectorSpace::Vector<T>& b, int* ipiv );
template <class T> inline void posv( VectorSpace::Matrix<T>& A, VectorSpace::Vector<T>& b );

namespace VectorSpace {

template <class T> class Vector
{
  private:
    T* _data;
    unsigned _N;
    bool _ref;

  public:
    typedef T* iterator;
    typedef const T* const_iterator;

    Vector();
    Vector( unsigned n, bool skipZeroInit = false );
    Vector( unsigned n, T init );
    Vector( unsigned n, T* ptr );
    Vector( const Vector& other );
    Vector( std::initializer_list<T> initList );
    ~Vector();
    void operator=( const Vector& other );
    void operator=( const std::vector<T>& other );

    T& operator[]( unsigned i );
    const T& operator[]( unsigned i ) const;
    T& operator[]( int i );
    const T& operator[]( int i ) const;

    Vector operator+( const Vector& other ) const;
    Vector operator-( const Vector& other ) const;
    Vector operator*( const Vector& other ) const;
    Vector operator/( const Vector& other ) const;

    Vector operator+( const T& scalar ) const;
    Vector operator-( const T& scalar ) const;
    Vector operator*( const T& scalar ) const;
    Vector operator/( const T& scalar ) const;

    void operator+=( const Vector& other );
    void operator-=( const Vector& other );
    void operator*=( const Vector& other );
    void operator/=( const Vector& other );

    void operator+=( const T& scalar );
    void operator-=( const T& scalar );
    void operator*=( const T& scalar );
    void operator/=( const T& scalar );

    unsigned size() const;
    void reset();
    void resize( unsigned newSize );
    iterator begin();
    iterator end();

    T inner( const Vector<T>& a ) const;

    friend void gesv<>( const Matrix<T>& A, Vector<T>& b, int* ipiv );
    friend void posv<>( Matrix<T>& A, Vector<T>& b );
};

template <class T> Vector<T>::Vector() : _data( nullptr ), _N( 0 ), _ref( false ) {}

template <class T> Vector<T>::Vector( unsigned n, bool skipZeroInit ) : _N( n ), _ref( false ), _data( nullptr ) {
    //_data = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    _data = new T[_N];
    if( !skipZeroInit ) {
        for( unsigned i = 0; i < _N; ++i ) {
            _data[i] = 0.0;
        }
    }
}

template <class T> Vector<T>::Vector( unsigned n, T init ) : _N( n ), _ref( false ), _data( nullptr ) {
    //_data = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = init;
    }
}

template <class T> Vector<T>::Vector( unsigned n, T* ptr ) : _N( n ), _ref( true ) { _data = ptr; }

template <class T> Vector<T>::~Vector() {
    if( !_ref && _data ) {
        // free( _data );
        delete[] _data;
    }
}

template <class T> Vector<T>::Vector( const Vector& other ) : _N( other._N ), _ref( false ), _data( nullptr ) {
    //_data = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = other._data[i];
    }
}

template <class T> Vector<T>::Vector( std::initializer_list<T> initList ) : _N( initList.size() ), _ref( false ), _data( nullptr ) {
    //_data        = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    _data        = new T[_N];
    auto listPtr = initList.begin();
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = *listPtr;
        ++listPtr;
    }
}

template <class T> void Vector<T>::operator=( const Vector& other ) {
    if( _data == nullptr ) {
        _N    = other._N;
        _ref  = false;
        _data = new T[_N];
        //_data = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    }
    assert( _N == other._N );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] = other._data[i];
    }
}

template <class T> void Vector<T>::operator=( const std::vector<T>& other ) {
    if( _data == nullptr ) {
        _N    = other.size();
        _ref  = false;
        _data = new T[_N];
        //_data = static_cast<T*>( malloc( _N * sizeof( T ) ) );
    }
    assert( _N == other.size() );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] = other[i];
    }
}

template <class T> T& Vector<T>::operator[]( unsigned i ) {
    assert( i < _N );
    return _data[i];
}

template <class T> const T& Vector<T>::operator[]( unsigned i ) const {
    assert( i < _N );
    return _data[i];
}

template <class T> T& Vector<T>::operator[]( int i ) {
    assert( i < _N );
    return _data[i];
}

template <class T> const T& Vector<T>::operator[]( int i ) const {
    assert( i < _N );
    return _data[i];
}

template <class T> Vector<T> Vector<T>::operator+( const Vector& other ) const {
    Vector<T> res( _N, true );
    assert( other._N == _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] + other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator-( const Vector& other ) const {
    assert( other._N == _N );
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] - other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator*( const Vector& other ) const {
    assert( other._N == _N );
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] * other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator/( const Vector& other ) const {
    assert( other._N == _N );
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] / other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator+( const T& scalar ) const {
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] + scalar;
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator-( const T& scalar ) const {
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] - scalar;
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator*( const T& scalar ) const {
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] * scalar;
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator/( const T& scalar ) const {
    Vector<T> res( _N, true );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] / scalar;
    }
    return res;
}

template <class T> void Vector<T>::operator+=( const Vector& other ) {
    assert( other._N == _N );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] += other._data[i];
    }
}

template <class T> void Vector<T>::operator-=( const Vector& other ) {
    assert( other._N == _N );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] -= other._data[i];
    }
}

template <class T> void Vector<T>::operator*=( const Vector& other ) {
    assert( other._N == _N );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] *= other._data[i];
    }
}

template <class T> void Vector<T>::operator/=( const Vector& other ) {
    assert( other._N == _N );
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] /= other._data[i];
    }
}

template <class T> void Vector<T>::operator+=( const T& scalar ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] += scalar;
    }
}

template <class T> void Vector<T>::operator-=( const T& scalar ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] -= scalar;
    }
}

template <class T> void Vector<T>::operator*=( const T& scalar ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] *= scalar;
    }
}

template <class T> void Vector<T>::operator/=( const T& scalar ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] /= scalar;
    }
}

template <class T> unsigned Vector<T>::size() const { return _N; }

template <class T> void Vector<T>::reset() {
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = 0.0;
    }
}

template <class T> void Vector<T>::resize( unsigned newSize ) {
    if( _data != nullptr ) {
        T* _newData = new T[newSize];
        // T* _newData = static_cast<T*>( malloc( newSize * sizeof( T ) ) );
        for( unsigned i = 0; i < std::min( _N, newSize ); ++i ) {
            _newData[i] = _data[i];
        }
        delete[] _data;
        // free( _data );
        _N    = newSize;
        _data = _newData;
    }
    else {
        _N    = newSize;
        _data = new T[newSize];
        //_data = static_cast<T*>( malloc( newSize * sizeof( T ) ) );
    }
}

template <class T> T* Vector<T>::begin() { return &_data[0]; }

template <class T> T* Vector<T>::end() { return &_data[_N]; }

template <class T> T Vector<T>::inner( const Vector<T>& a ) const {
    double tmp = 0.0;
    assert( a.size() == _N );
    for( unsigned i = 0; i < a.size(); ++i ) {
        tmp += a[i] * _data[i];
    }
    return tmp;
}

}    // namespace VectorSpace
