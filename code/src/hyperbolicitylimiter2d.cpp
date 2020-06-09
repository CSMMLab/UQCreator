#include "hyperbolicitylimiter2d.h"
#include "mathtools.h"

HyperbolicityLimiter2D::HyperbolicityLimiter2D( Settings* settings )
    : Closure( settings ), _gamma( _settings->GetGamma() ), _lambda( _settings->GetFilterStrength() ) {
    _alpha = 1.0;
    if( _lambda < 0 ) _lambda = 0.0;
    unsigned maxDegree = _settings->GetMaxDegree();
    _filterFunction    = Vector( _settings->GetNTotal(), 1.0 );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index =
                    unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
                _filterFunction[i] *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) );
            }
        }
    }
}

HyperbolicityLimiter2D::~HyperbolicityLimiter2D() {}

void HyperbolicityLimiter2D::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void HyperbolicityLimiter2D::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix HyperbolicityLimiter2D::U( const Matrix& Lambda ) { return Lambda; }

void HyperbolicityLimiter2D::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void HyperbolicityLimiter2D::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector rho( _nQTotalForRef[refLevel], 0.0 );
    Vector m1( _nQTotalForRef[refLevel], 0.0 );
    Vector m2( _nQTotalForRef[refLevel], 0.0 );
    Vector E( _nQTotalForRef[refLevel], 0.0 );
    Vector p2( _nQTotalForRef[refLevel], 0.0 );
    Vector t1( _nQTotalForRef[refLevel], 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    // save zero order moments (realizable solution)
    double rhoTilde = u( 0, 0 );
    double m1Tilde  = u( 1, 0 );
    double m2Tilde  = u( 2, 0 );
    double ETilde   = u( 3, 0 );

    // save zero order moment as full moment vector for output
    Matrix u2( _nStates, nTotal );
    u2( 0, 0 ) = u( 0, 0 );
    u2( 1, 0 ) = u( 1, 0 );
    u2( 2, 0 ) = u( 2, 0 );
    u2( 3, 0 ) = u( 3, 0 );

    double t1Max = -1000;

    double p1 = ( _gamma - 1.0 ) * ( ETilde - 0.5 * pow( m1Tilde, 2 ) / rhoTilde - 0.5 * pow( m2Tilde, 2 ) / rhoTilde );

    // Question: Is e a scalar in our case?
    double e = MathTools::min( 1e-10, MathTools::min( rhoTilde, p1 ) );

    // determine t1Max
    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        Vector uKinetic = EvaluateLambda( u, k, nTotal );

        rho[k] = uKinetic[0];
        m1[k]  = uKinetic[1];
        m2[k]  = uKinetic[2];
        E[k]   = uKinetic[3];

        t1[k] = ( e - rho[k] ) / ( rhoTilde - rho[k] );
        if( t1[k] > 1.0 || t1[k] < 0.0 || !std::isfinite( t1[k] ) ) t1[k] = 0.0;

        if( t1[k] > t1Max ) t1Max = t1[k];
    }

    double t2Max = -1000;

    // determine t2Max
    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        // rho[k] = t1Max * rhoTilde + ( 1.0 - t1Max ) * rho[k];    // why do we add this (taken from scalar code)

        // a* t ^ 2 + b* t + c = 0
        double a = ( E[k] - ETilde ) * ( rho[k] - rhoTilde ) - 0.5 * pow( m1[k] - m1Tilde, 2 ) - 0.5 * pow( m2[k] - m2Tilde, 2 );
        double b = m1[k] * ( m1[k] - m1Tilde ) + m2[k] * ( m2[k] - m2Tilde ) - rho[k] * ( E[k] - ETilde ) - E[k] * ( rho[k] - rhoTilde );
        double c = E[k] * rho[k] - 0.5 * pow( m1[k], 2 ) - 0.5 * pow( m2[k], 2 ) - e;

        // Citardauq Formula as stable p-q formula
        double q = -0.5 * ( b + MathTools::sign( b ) * sqrt( pow( b, 2 ) - 4.0 * a * c ) );
        if( pow( b, 2 ) - 4.0 * a * c < 0 ) std::cerr << "[RealizabilityLimiter]: Imaginary unit detected!" << std::endl;
        double t2a = q / a;
        double t2b = c / q;

        if( t2a > 1.0 || t2a < 0.0 || !std::isfinite( t2a ) ) t2a = 0.0;
        if( t2b > 1.0 || t2b < 0.0 || !std::isfinite( t2b ) ) t2b = 0.0;

        double t2MaxTmp = MathTools::max( t2a, t2b );
        if( t2MaxTmp > t2Max ) t2Max = t2MaxTmp;
    }

    double theta = MathTools::max( t2Max, t1Max );

    lambda = theta * u2 + ( 1.0 - theta ) * u;
}

// taken from 1D limiter and extended to 2D
/*
void HyperbolicityLimiter2D::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector rho( _nQTotalForRef[refLevel], 0.0 );
    Vector m1( _nQTotalForRef[refLevel], 0.0 );
    Vector m2( _nQTotalForRef[refLevel], 0.0 );
    Vector E( _nQTotalForRef[refLevel], 0.0 );
    Vector p( _nQTotalForRef[refLevel], 0.0 );
    Vector t1( _nQTotalForRef[refLevel], 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    // save zero order moments (realizable solution)
    double rhoTilde = u( 0, 0 );
    double m1Tilde  = u( 1, 0 );
    double m2Tilde  = u( 2, 0 );
    double ETilde   = u( 3, 0 );

    // save zero order moment as full moment vector for output
    Matrix u2( _nStates, nTotal );
    u2( 0, 0 ) = u( 0, 0 );
    u2( 1, 0 ) = u( 1, 0 );
    u2( 2, 0 ) = u( 2, 0 );
    u2( 3, 0 ) = u( 3, 0 );

    double t1Max = -1000;

    double pTilde = ( _gamma - 1.0 ) * ( ETilde - 0.5 * pow( m1Tilde, 2 ) / rhoTilde - 0.5 * pow( m2Tilde, 2 ) / rhoTilde );

    // Question: Is e a scalar in our case?
    double e = MathTools::min( 1e-10, MathTools::min( rhoTilde, pTilde ) );

    // determine t1Max
    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        Vector uKinetic = EvaluateLambda( u, k, nTotal );

        rho[k] = uKinetic[0];
        m1[k]  = uKinetic[1];
        m2[k]  = uKinetic[2];
        E[k]   = uKinetic[2];

        t1[k] = ( rho[k] - e ) / ( rho[k] - rhoTilde );
        if( t1[k] > 1.0 || t1[k] < 0.0 || !std::isfinite( t1[k] ) ) t1[k] = 0.0;

        if( t1[k] > t1Max ) t1Max = t1[k];
    }

    double t2Max = -1000;

    // determine t2Max
    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        rho[k] = t1Max * rhoTilde + ( 1.0 - t1Max ) * rho[k];    // why should we do this??

        // a* t ^ 2 + b* t + c = 0
        double a = ( ETilde - E[k] ) * ( rhoTilde - rho[k] ) - 0.5 * pow( m1[k] - m1Tilde, 2 ) - 0.5 * pow( m2[k] - m2Tilde, 2 );
        double b = m1[k] * ( m1[k] - m1Tilde ) - E[k] * ( rho[k] - rhoTilde ) - rho[k] * ( E[k] - ETilde ) + m2[k] * ( m2[k] - m2Tilde );
        double c = E[k] * rho[k] - 0.5 * pow( m1[k], 2 ) - 0.5 * pow( m2[k], 2 ) - e;

        double q = -0.5 * ( b + MathTools::sign( b ) * sqrt( pow( b, 2 ) - 4.0 * a * c ) );
        if( pow( b, 2 ) - 4.0 * a * c < 0 ) std::cerr << "[RealizabilityLimiter]: Imaginary unit detected!" << std::endl;
        // how can I check this in C++ ?
        // q        = real( q );
        double t2a = q / a;
        double t2b = c / q;

        if( t2a > 1.0 || t2a < 0.0 || !std::isfinite( t2a ) ) t2a = 0.0;
        if( t2b > 1.0 || t2b < 0.0 || !std::isfinite( t2b ) ) t2b = 0.0;

        double t2MaxTmp = MathTools::max( t2a, t2b );
        if( t2MaxTmp > t2Max ) t2Max = t2MaxTmp;
    }
    double theta = t2Max;
    if( t1Max > t2Max ) theta = t1Max;
    lambda = theta * u2 + ( 1 - theta ) * u;
}*/

void HyperbolicityLimiter2D::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) { this->SolveClosure( lambda, u, refLevel ); }
