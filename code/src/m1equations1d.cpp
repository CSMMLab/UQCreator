#include "m1equations1d.h"
#include "mathtools.h"

M1Equations1D::M1Equations1D( Settings* settings ) : Problem( settings ) {
    _nStates = 2;
    _settings->SetNStates( _nStates );
    _settings->SetSource( true );
    _sigmaA = 0.0;    // absorption coefficient
    _sigmaS = 1.0;    // scattering coefficient
    _sigmaT = _sigmaA + _sigmaS;
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[M1Equations1D] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}
M1Equations1D::~M1Equations1D() {}

Vector M1Equations1D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    return 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
}

Matrix M1Equations1D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix M1Equations1D::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );

    double alpha = Bisection( -500.0, 100.1, u[1] / u[0] );    // ComputeAlpha( u[1] / u[0] );
    double u2Du0 = 4.0 / ( alpha - alpha * exp( 2.0 * alpha ) ) - 2.0 / alpha + 2.0 / std::pow( alpha, 2 ) + 1.0;

    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u2Du0 * u[0];

    return flux;
}

double M1Equations1D::ComputeAlpha( const double u1Du0 ) const {
    double H;
    double alpha   = -1e-5;
    double epsilon = 1e-3;

    double g = MathTools::coth( alpha ) - 1.0 / alpha - u1Du0;

    while( fabs( g ) > epsilon ) {
        H     = 1.0 / pow( alpha, 2 ) - MathTools::csch( alpha );
        alpha = alpha - g / H;
        g     = MathTools::coth( alpha ) - 1.0 / alpha - u1Du0;
        // std::cout << "alpha = " << alpha << " ; res = " << fabs( g ) << std::endl;
    }
    return alpha;
}

double M1Equations1D::RootFun( const double alpha, const double u1Du0 ) const { return MathTools::coth( alpha ) - 1.0 / alpha - u1Du0; }

double M1Equations1D::Bisection( double alphaA, double alphaB, const double u1Du0 ) const {
    if( RootFun( alphaA, u1Du0 ) * RootFun( alphaB, u1Du0 ) >= 0 ) {
        std::cerr << "[Bisection] Incorrect a and b with vals " << RootFun( alphaA, u1Du0 ) << ", " << RootFun( alphaB, u1Du0 ) << std::endl;
        exit( EXIT_FAILURE );
        return -1.0;
    }

    double alphaC = alphaA;
    double e      = 1e-10;

    while( std::abs( alphaB - alphaA ) >= e ) {
        alphaC = ( alphaA + alphaB ) / 2;
        if( !std::isfinite( RootFun( alphaC, u1Du0 ) ) ) {
            std::cerr << "[Bisection] Infinite RootFunction" << std::endl;
            exit( EXIT_FAILURE );
        }
        else if( RootFun( alphaC, u1Du0 ) == 0.0 ) {
            return alphaC;
        }
        else if( RootFun( alphaC, u1Du0 ) * RootFun( alphaA, u1Du0 ) < 0 ) {
            alphaB = alphaC;
        }
        else {
            alphaA = alphaC;
        }
    }
    return alphaC;
}

Matrix M1Equations1D::F( const Matrix& u ) {
    _log->error( "[euler2d] Flux not implemented" );
    exit( EXIT_FAILURE );
    return 0.5 * pow( u, 2 );
}

Matrix M1Equations1D::Source( const Matrix& uQ ) const {
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = static_cast<unsigned>( uQ.columns() );
    Vector g( nStates, 0.0 );
    g[0] = 1.0;    // 2 * M_PI;
    Matrix y( nStates, Nq, 0.0 );
    for( unsigned s = 0; s < nStates; ++s ) {
        for( unsigned k = 0; k < Nq; ++k ) {
            y( s, k ) = -_sigmaA * uQ( s, k ) - _sigmaS * ( 1.0 - g[s] ) * uQ( s, k );
            // _sigmaA* uQ( s, k ) + _sigmaS * (1 - g[s]) * uQ( s, k );
        }
    }
    return y;
}

double M1Equations1D::ComputeDt( const Matrix& u, double dx, unsigned level ) const {

    double cfl = _settings->GetCFL();

    double maxVelocity = 1.0;

    return ( cfl * dx ) / maxVelocity;
}

Vector M1Equations1D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double x0      = 0.0;
    double s2      = 1.0 * std::pow( 0.01, 2 );    // std::pow( 0.03, 2 );
    double floor   = 1e-4;
    _sigma         = _settings->GetSigma();
    double sigmaXi = _sigma[0] * xi[0];

    y[0] = std::fmax(
        floor, pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) * exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );

    return y;
}

Vector M1Equations1D::LoadIC( const Vector& x, const Vector& xi ) { return Vector( 1 ); }
