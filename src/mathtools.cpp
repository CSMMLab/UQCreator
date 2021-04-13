#include "mathtools.h"

double MathTools::Pythag( const double a, const double b ) {
    double absa = std::fabs( a ), absb = std::fabs( b );
    return ( absa > absb ? absa * std::sqrt( 1.0 + std::pow( absb / absa, 2 ) )
                         : ( absb == 0.0 ? 0.0 : absb * std::sqrt( 1.0 + std::pow( absa / absb, 2 ) ) ) );
}

unsigned MathTools::BinomialCoefficient( unsigned n, unsigned k ) {
    if( k > n ) return 0;

    if( n - k < k ) k = n - k;

    unsigned r = 1;

    for( unsigned d = 1; d <= k; d++ ) {
        unsigned mult = n;

        bool divided = true;

        if( mult % d == 0 )
            mult /= d;
        else if( r % d == 0 )
            r /= d;
        else
            divided = false;

        const unsigned r_mult = r * mult;
        if( r_mult / mult != r ) throw std::overflow_error( "Overflow" );

        r = r_mult;

        if( !divided ) r /= d;

        n--;
    }

    return r;
}

std::pair<Vector, Matrix> MathTools::ComputeEigenValTriDiagMatrix( const Matrix& mat ) {
    unsigned n = mat.rows();

    Vector d( n, 0.0 ), e( n, 0.0 );
    Matrix z( n, n, 0.0 );
    for( unsigned i = 0; i < n; ++i ) {
        d[i]          = mat( i, i );
        z( i, i )     = 1.0;
        i == 0 ? e[i] = 0.0 : e[i] = mat( i, i - 1 );
    }

    int m, l, iter, i, k;
    m = l = iter = i = k = 0;
    double s, r, p, g, f, dd, c, b;
    s = r = p = g = f = dd = c = b = 0.0;
    const double eps               = std::numeric_limits<double>::epsilon();
    for( i = 1; i < static_cast<int>( n ); i++ ) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    for( l = 0; l < static_cast<int>( n ); l++ ) {
        iter = 0;
        do {
            for( m = l; m < static_cast<int>( n ) - 1; m++ ) {
                dd = std::fabs( d[m] ) + std::fabs( d[m + 1] );
                if( std::fabs( e[m] ) <= eps * dd ) break;
            }
            if( m != l ) {
                if( iter++ == 30 ) throw( "[solveTriDiagEWP]: Too many iterations" );
                g = ( d[l + 1] - d[l] ) / ( 2.0 * e[l] );
                r = Pythag( g, 1.0 );
                g = d[m] - d[l] + e[l] / ( g + std::copysign( r, g ) );
                s = c = 1.0;
                p     = 0.0;
                for( i = m - 1; i >= l; i-- ) {
                    f        = s * e[i];
                    b        = c * e[i];
                    e[i + 1] = ( r = Pythag( f, g ) );
                    if( r == 0.0 ) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s        = f / r;
                    c        = g / r;
                    g        = d[i + 1] - p;
                    r        = ( d[i] - g ) * s + 2.0 * c * b;
                    d[i + 1] = g + ( p = s * r );
                    g        = c * r - b;
                    for( k = 0; k < static_cast<int>( n ); k++ ) {
                        f = z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) + 1 );
                        z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) + 1 ) =
                            s * z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) + c * f;
                        z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) =
                            c * z( static_cast<unsigned>( k ), static_cast<unsigned>( i ) ) - s * f;
                    }
                }
                if( r == 0.0 && i >= l ) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while( m != l );
    }
    return std::make_pair( d, z );
}

int MathTools::Factorial( int i ) {
    int fact = 1;
    for( int n = 1; n <= i; n++ ) {
        fact *= i;
    }
    return fact;
}

double MathTools::max( double a, double b ) {
    if( a > b )
        return a;
    else
        return b;
}

double MathTools::min( double a, double b ) {
    if( b < a )
        return b;
    else
        return a;
}

double MathTools::sign( double a ) {
    if( a > 0 )
        return 1;
    else
        return -1;
}

double MathTools::csch( const double x ) {
    // return 1.0 / sinh( x );
    return fabs( x ) / ( x * fabs( sinh( x ) ) );
}

double MathTools::coth( const double x ) { return sinh( 2.0 * x ) / ( cosh( 2.0 * x ) - 1.0 ); }
