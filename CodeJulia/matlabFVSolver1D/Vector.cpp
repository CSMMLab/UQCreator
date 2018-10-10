#include "Vector.h"

Vector::Vector( const unsigned N ) : _n( N ) { _a = new double[_n]; }

Vector::~Vector() { delete[] _a; }

Vector::Vector( const Vector& R ) : _n( R._n ) {
    if( &R == this ) return;
    _a = new double[_n];
    for( unsigned i = 0; i < R._n; ++i ) {
        _a[i] = R._a[i];
    }
    return;
}

Vector& Vector::operator=( const Vector& R ) {
    if( &R == this ) return *this;
    for( unsigned i = 0; i < _n; ++i ) {
        _a[i] = R._a[i];
    }
    return *this;
}

double& Vector::operator[]( const unsigned i ) {
    assert( i < _n );
    return _a[i];
}

const double& Vector::operator[]( const unsigned i ) const {
    assert( i < _n );
    return _a[i];
}

Vector Vector::operator*( const double Scalar ) const {
    Vector r( _n );
    for( unsigned i = 0; i < _n; ++i ) {
        r._a[i] = _a[i] * Scalar;
    }
    return r;
}

Vector Vector::operator+( const Vector& B ) const {
    assert( _n == B._n );
    Vector r( _n );
    for( unsigned i = 0; i < _n; ++i ) {
        r._a[i] = B._a[i] + _a[i];
    }
    return r;
}

Vector Vector::operator-( const Vector& B ) const {
    assert( _n == B._n );
    Vector r( _n );
    for( unsigned i = 0; i < _n; ++i ) {
        r._a[i] = _a[i] - B._a[i];
    }
    return r;
}

void Vector::Show() const {
    std::cout << std::endl;
    for( unsigned i = 0; i < _n; ++i ) {
        std::cout << _a[i] << " ";
    }
    std::cout << std::endl;
    return;
}

void Vector::Add( Vector* a ) {
    assert( _n == a->Dimension() );
    for( unsigned i = 0; i < _n; ++i ) {
        _a[i] += ( *a )[i];
    }
}

void Vector::Subtract( Vector* a ) {
    assert( _n == a->Dimension() );
    for( unsigned i = 0; i < _n; ++i ) {
        _a[i] -= ( *a )[i];
    }
}

void Vector::Multiply( const double Scalar ) {
    for( unsigned i = 0; i < _n; ++i ) {
        _a[i] *= Scalar;
    }
}

void Vector::MakeZero() {
    for( unsigned i = 0; i < _n; ++i ) _a[i] = 0;
}
