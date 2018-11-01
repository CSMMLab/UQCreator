#include <cmath>
#include <initializer_list>
#include <iostream>
#include <sstream>

namespace VectorSpace {
template <class T> class Vector;
template <class T> class Matrix;
}    // namespace VectorSpace

template <class T> inline void gesv( VectorSpace::Matrix<T>& A, VectorSpace::Vector<T>& b, int* ipiv );

namespace VectorSpace {

template <class T> class Vector
{
  private:
    unsigned _N;
    T* _data;
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

    T& operator[]( unsigned i );
    const T& operator[]( unsigned i ) const;
    Vector operator+( const Vector& other ) const;
    Vector operator-( const Vector& other ) const;
    Vector operator*( const Vector& other ) const;
    Vector operator*( const T& scalar ) const;
    Vector operator/( const Vector& other ) const;
    void operator=( const Vector& other );
    void operator+=( const Vector& other );
    void operator-=( const Vector& other );
    void operator*=( const Vector& other );
    void operator*=( const T& scalar );
    void operator/=( const Vector& other );

    unsigned size() const;
    void reset();
    void resize( unsigned newSize );
    iterator begin();
    iterator end();

    friend void gesv<>( Matrix<T>& A, Vector<T>& b, int* ipiv );
};

template <class T> Vector<T>::Vector() : _data( nullptr ), _N( 0 ), _ref( false ) {}

template <class T> Vector<T>::Vector( unsigned n, bool skipZeroInit ) : _N( n ), _ref( false ) {
    _data = new T[_N];
    if( !skipZeroInit ) {
        for( unsigned i = 0; i < _N; ++i ) {
            _data[i] = 0.0;
        }
    }
}

template <class T> Vector<T>::Vector( unsigned n, T init ) : _N( n ), _ref( false ) {
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = init;
    }
}

template <class T> Vector<T>::Vector( unsigned n, T* ptr ) : _N( n ), _ref( true ) { _data = ptr; }

template <class T> Vector<T>::~Vector() {
    if( !_ref ) {
        delete[] _data;
    }
}

template <class T> Vector<T>::Vector( const Vector& other ) : _N( other._N ), _ref( false ) {
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = other._data[i];
    }
}

template <class T> Vector<T>::Vector( std::initializer_list<T> initList ) : _N( initList.size() ), _ref( false ) {
    _data        = new T[_N];
    auto listPtr = initList.begin();
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = *listPtr;
        ++listPtr;
    }
}

template <class T> T& Vector<T>::operator[]( unsigned i ) { return _data[i]; }

template <class T> const T& Vector<T>::operator[]( unsigned i ) const { return _data[i]; }

template <class T> Vector<T> Vector<T>::operator+( const Vector& other ) const {
    Vector<T> res( _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] + other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator-( const Vector& other ) const {
    Vector<T> res( _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] - other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator*( const Vector& other ) const {
    Vector<T> res( _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] * other._data[i];
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator*( const T& scalar ) const {
    Vector<T> res( _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] * scalar;
    }
    return res;
}

template <class T> Vector<T> Vector<T>::operator/( const Vector& other ) const {
    Vector<T> res( _N );
    for( unsigned i = 0; i < _N; ++i ) {
        res[i] = this->_data[i] / other._data[i];
    }
    return res;
}

template <class T> void Vector<T>::operator=( const Vector& other ) {
    if( _data == nullptr ) {
        _N    = other._N;
        _ref  = false;
        _data = new T[_N];
    }
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] = other._data[i];
    }
}

template <class T> void Vector<T>::operator+=( const Vector& other ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] += other._data[i];
    }
}

template <class T> void Vector<T>::operator-=( const Vector& other ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] -= other._data[i];
    }
}

template <class T> void Vector<T>::operator*=( const Vector& other ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] *= other._data[i];
    }
}

template <class T> void Vector<T>::operator*=( const T& scalar ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] *= scalar;
    }
}

template <class T> void Vector<T>::operator/=( const Vector& other ) {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] /= other._data[i];
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
        for( unsigned i = 0; i < std::min( _N, newSize ); ++i ) {
            _newData[i] = _data[i];
        }
        delete[] _data;
        _N    = newSize;
        _data = _newData;
    }
    else {
        _N    = newSize;
        _data = new T[newSize];
    }
}

template <class T> T* Vector<T>::begin() { return &_data[0]; }

template <class T> T* Vector<T>::end() { return &_data[_N]; }

}    // namespace VectorSpace
