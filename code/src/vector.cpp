#include "vector.h"

namespace VectorSpace {

template <class T> Vector<T>::Vector() : _N( 0 ) {}

template <class T> Vector<T>::Vector( unsigned n, bool skipZeroInit ) : _N( n ) {
    _data = new T[_N];
    if( !skipZeroInit ) {
        for( unsigned i = 0; i < _N; ++i ) {
            _data[i] = 0.0;
        }
    }
}

template <class T> Vector<T>::Vector( unsigned n, T init ) : _N( n ) {
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = init;
    }
}

template <class T> Vector<T>::~Vector() { delete[] _data; }

template <class T> Vector<T>::Vector( const Vector& other ) {
    _N    = other._N;
    _data = new T[_N];
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = other._data[i];
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

template <class T> Vector<T> Vector<T>::operator*( const double& scalar ) const {
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

template <class T> void Vector<T>::operator=( const Vector& other ) const {
    for( unsigned i = 0; i < _N; ++i ) {
        this->_data[i] = other._data[i];
    }
}

template <class T> unsigned Vector<T>::size() const { return _N; }

template <class T> void Vector<T>::reset() {
    for( unsigned i = 0; i < _N; ++i ) {
        _data[i] = 0.0;
    }
}

template <class T> void Vector<T>::resize( unsigned newSize ) {
    T* _newData = new T[newSize];
    for( unsigned i = 0; i < std::min( _N, newSize ); ++i ) {
        _newData[i] = _data[i];
    }
    delete[] _data;
    _data = _newData;
}

template <class T> T* Vector<T>::begin() { return &_data[0]; }

template <class T> T* Vector<T>::end() { return &_data[_N]; }

}    // namespace VectorSpace
