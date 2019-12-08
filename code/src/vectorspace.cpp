#include <iostream>
#include <spdlog/spdlog.h>

#include "matrix.cpp"

/* Complex datatype */
struct _fcomplex {
    float re, im;
};
typedef struct _fcomplex fcomplex;

extern "C" {
void dgesv_( int* n, int* nrhs, double* A, int* lda, int* ipiv, double* b, int* ldb, int* info );
void dposv_( char* uplo, int* n, int* nrhs, double* A, int* lda, double* b, int* ldb, int* info );
void cgeev_( char* jobvl,
             char* jobvr,
             int* n,
             fcomplex* a,
             int* lda,
             fcomplex* w,
             fcomplex* vl,
             int* ldvl,
             fcomplex* vr,
             int* ldvr,
             fcomplex* work,
             int* lwork,
             float* rwork,
             int* info );
}

// MATRIX //////////////////////////////////////////////////////////////////////////

template <class T> VectorSpace::Vector<T> column( VectorSpace::Matrix<T>& mat, unsigned i ) {
    return VectorSpace::Vector<T>( mat._rows, &mat._data[i * mat._rows] );
}

template <class T> VectorSpace::Vector<T> column( const VectorSpace::Matrix<T>& mat, unsigned i ) {
    return VectorSpace::Vector<T>( mat._rows, &mat._data[i * mat._rows] );
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
    VectorSpace::Matrix<T> res( mat.columns(), mat.rows() );
    for( unsigned j = 0; j < mat.rows(); ++j ) {
        for( unsigned i = 0; i < mat.columns(); ++i ) {
            res( i, j ) = mat( j, i );
        }
    }
    return res;
}

template <class T> std::ostream& operator<<( std::ostream& os, const VectorSpace::Matrix<T>& a ) {
    for( unsigned i = 0; i < a.rows(); ++i ) {
        os << "(\t";
        for( unsigned j = 0; j < a.columns(); ++j ) {
            os << a( i, j ) << "\t";
        }
        os << ")" << std::endl;
    }
    return os;
}

// VECTOR //////////////////////////////////////////////////////////////////////////

template <class T> double dot( const VectorSpace::Vector<T>& a, const VectorSpace::Vector<T>& b ) {
    double res = 0.0;
    for( unsigned i = 0; i < a.size(); ++i ) {
        res += a[i] * b[i];
    }
    return res;
}

template <class T> VectorSpace::Matrix<T> outer( const VectorSpace::Vector<T>& a, const VectorSpace::Vector<T>& b ) {
    VectorSpace::Matrix<T> res( a.size(), a.size() );
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
    for( unsigned i = 0; i < a.size(); ++i ) {
        os << "( " << a[i] << " )" << std::endl;
    }
    return os;
}

template <class T> double l2Norm( const VectorSpace::Vector<T>& a ) {
    double res = 0;
    for( unsigned i = 0; i < a.size(); ++i ) {
        double tmp = static_cast<double>( a[i] );
        res += std::fabs( tmp ) * std::fabs( tmp );
    }
    return std::sqrt( res );
}

template <class T> double norm( const VectorSpace::Vector<T>& a ) { return l2Norm<T>( a ); }

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

// SOLVER //////////////////////////////////////////////////////////////////////////

template <class T> inline void gesv( const VectorSpace::Matrix<T>& A, VectorSpace::Vector<T>& b, int* ipiv ) {
    int n( A.rows() );
    int nrhs( 1 );
    int lda( A.columns() );
    int ldb( b.size() );
    int info( 0 );

    if( n == 0 ) {
        return;
    }

    VectorSpace::Matrix<T> ATmp = A;

    dgesv_( &n, &nrhs, ATmp._data, &lda, ipiv, b._data, &ldb, &info );

    if( info > 0 ) {
        auto log = spdlog::get( "event" );
        log->error( "[gesv] Inversion of singular matrix failed" );
        exit( EXIT_FAILURE );
    }
}

template <class T> inline void posv( VectorSpace::Matrix<T>& A, VectorSpace::Vector<T>& b ) {
    int n( A.rows() );
    int nrhs( 1 );
    int lda( A.columns() );
    int ldb( b.size() );
    int info( 0 );
    char lu( 'U' );

    if( n == 0 ) {
        return;
    }

    dposv_( &lu, &n, &nrhs, A._data, &lda, b._data, &ldb, &info );

    if( info > 0 ) {
        auto log = spdlog::get( "event" );
        log->error( "[posv] Inversion of singular matrix failed" );
        exit( EXIT_FAILURE );
    }
}

template <class T>
inline void cgeev( const VectorSpace::Matrix<T>& A, VectorSpace::Matrix<T>& VL, VectorSpace::Matrix<T>& VR, VectorSpace::Matrix<T>& W ) {

    /* Locals */
    int N = A.columns();
    int n = A.columns(), lda = A.columns(), ldvl = A.columns(), ldvr = A.columns(), info, lwork;
    fcomplex wkopt;
    fcomplex* work;
    /* Local arrays */
    /* rwork dimension should be at least 2*n */
    float rwork[2 * N];
    fcomplex w[N], vl[N * N], vr[N * N];
    fcomplex a[N * N];
    for( unsigned i = 0; i < N; ++i ) {
        for( unsigned j = 0; j < N; ++j ) {
            a[i + j * N] = {A( i, j ), 0.0};
        }
    }

    /* Executable statements */
    printf( " CGEEV Example Program Results\n" );
    /* Query and allocate the optimal workspace */
    lwork = -1;
    cgeev_( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, rwork, &info );
    lwork = (int)wkopt.re;
    work  = (fcomplex*)malloc( lwork * sizeof( fcomplex ) );
    /* Solve eigenproblem */
    cgeev_( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info );
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    /* Print eigenvalues */
    // print_matrix( "Eigenvalues", 1, n, w, 1 );
    /* Print left eigenvectors */
    // print_matrix( "Left eigenvectors", n, n, vl, ldvl );
    /* Print right eigenvectors */
    // print_matrix( "Right eigenvectors", n, n, vr, ldvr );
    for( unsigned i = 0; i < N; ++i ) {
        for( unsigned j = 0; j < N; ++j ) {
            VL( i, j ) = vl[i + j * N].re;
            VR( i, j ) = vr[i + j * N].re;
        }
        W( i, i ) = w[i].re;
    }
}

// multiplies PE part of A and x and saves result on b
template <class T>
void multOnPENoReset(
    const VectorSpace::Matrix<T>& X, const VectorSpace::Matrix<T>& A, VectorSpace::Matrix<T>& B, unsigned kStart, unsigned kEnd, unsigned nTotal ) {
    // B.reset();
    for( unsigned s = 0; s < B.rows(); ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned k = kStart; k <= kEnd; ++k ) {
                B( s, i ) += X( s, k - kStart ) * A( k, i );    // TODO do not choose PhiTiledeTrans but PhiTiled s.t. A( j,k ) can be used
            }
        }
    }
}
