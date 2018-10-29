#include "matrix.h"

// MATRIX //////////////////////////////////////////////////////////////////////////

template <class T> VectorSpace::Vector<T&> row( VectorSpace::Matrix<T>& mat, unsigned i ) {
    VectorSpace::Vector<T&> res( mat.columns() );
    for( unsigned j = 0; j < mat.columns(); ++j ) {
        res[j] = mat( i, j );
    }
    return res;
}

template <class T> VectorSpace::Vector<T&> column( VectorSpace::Matrix<T>& mat, unsigned i ) {
    VectorSpace::Vector<T&> res( mat.rows() );
    for( unsigned j = 0; j < mat.rows(); ++j ) {
        res[j] = mat( j, i );
    }
    return res;
}

template <class T> VectorSpace::Vector<T> row( const VectorSpace::Matrix<T>& mat, unsigned i ) {
    VectorSpace::Vector<T> res( mat.columns() );
    for( unsigned j = 0; j < mat._columns; ++j ) {
        res[j] = mat( i, j );
    }
    return res;
}

template <class T> VectorSpace::Vector<T> column( const VectorSpace::Matrix<T>& mat, unsigned i ) {
    VectorSpace::Vector<T> res( mat.rows() );
    for( unsigned j = 0; j < mat.rows(); ++j ) {
        res[j] = mat( j, i );
    }
    return res;
}

template <class T> VectorSpace::Matrix<T> pow( const VectorSpace::Matrix<T>& a, unsigned pot ) {
    if( pot == 0 ) {
        return VectorSpace::Matrix<T>( a.rows(), a.columns(), 1.0 );
    }
    VectorSpace::Matrix<T> res( a );
    for( unsigned j = 0; j < pot - 1; ++j ) {
        res = res * a;
    }
    return res;
}

template <class T, class U> VectorSpace::Matrix<T> operator*( const U& scalar, const VectorSpace::Matrix<T>& mat ) {
    VectorSpace::Matrix<T> res( mat.rows(), mat.columns(), true );
    for( unsigned j = 0; j < mat.columns(); ++j ) {
        for( unsigned i = 0; i < mat.rows(); ++i ) {
            res( i, j ) = mat( i, j ) * scalar;
        }
    }
    return res;
}

template <class T> VectorSpace::Matrix<T> trans( const VectorSpace::Matrix<T>& mat ) {
    VectorSpace::Matrix<T> res( mat.rows(), mat.columns() );
    for( unsigned j = 0; j < mat.columns(); ++j ) {
        for( unsigned i = 0; i < mat.rows(); ++i ) {
            res( i, j ) = mat( j, i );
        }
    }
    return res;
}

// VECTOR //////////////////////////////////////////////////////////////////////////

template <class T> double dot( const VectorSpace::Vector<T>& a, const VectorSpace::Vector<T>& b ) {
    double res = 0.0;
    for( unsigned i = 0; i < a._N; ++i ) {
        res += a._data * b._data[i];
    }
    return res;
}

template <class T> VectorSpace::Matrix<T> outer( const VectorSpace::Vector<T>& a, const VectorSpace::Vector<T>& b ) {
    VectorSpace::Matrix<T> res( a.size(), a.size(), true );
    for( unsigned i = 0; i < a.size(); ++i ) {
        for( unsigned j = 0; j < b.size(); ++j ) {
            res( i, j ) += a[i] * b[j];
        }
    }
    return res;
}

template <class T> VectorSpace::Vector<T> pow( const VectorSpace::Vector<T>& a, unsigned pot ) {
    if( pot == 0 ) {
        return VectorSpace::Vector<T>( a._N, 1.0 );
    }

    VectorSpace::Vector<T> res( a );
    for( unsigned i = 0; i < a._N; ++i ) {
        for( unsigned j = 0; j < pot - 1; ++j ) {
            res[i] *= a._data[i];
        }
    }
}

template <class T> std::ostream& operator<<( std::ostream& os, const VectorSpace::Vector<T>& a ) {
    for( unsigned i = 0; i < a._N; ++i ) {
        os << "( " << a[i] << " )" << std::endl;
    }
    os << std::endl;
    return os;
}

template <class T> double l2Norm( const VectorSpace::Vector<T>& a ) {
    double res = 0;
    for( unsigned i = 0; i < a.size(); ++i ) {
        res += std::fabs( a[i] ) * std::fabs( a[i] );
    }
    return std::sqrt( res );
}

template <class T> double norm( const VectorSpace::Vector<T>& a ) { return l2Norm( a ); }

template <class T> VectorSpace::Vector<T> operator*( const double& scalar, const VectorSpace::Vector<T>& vec ) {
    VectorSpace::Vector<T> res( vec.size(), true );
    for( unsigned i = 0; i < vec.size(); ++i ) {
        res[i] = scalar * vec[i];
    }
    return res;
}

template <class T> VectorSpace::Vector<T> operator-( const VectorSpace::Vector<T>& vec ) {
    VectorSpace::Vector<T> res( vec.size(), true );
    for( unsigned i = 0; i < vec.size(); ++i ) {
        res[i] = -vec[i];
    }
    return res;
}
