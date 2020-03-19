#include "hyperbolicitylimiter2d.h"
#include "mathtools.h"

HyperbolicityLimiter2D::HyperbolicityLimiter2D( Settings* settings )
    : Closure( settings ), _gamma( _settings->GetGamma() ), _lambda( _settings->GetFilterStrength() ) {
    _alpha = 1.0;
    if( _lambda < 0 ) _lambda = 0.0;
    unsigned nMoments = _settings->GetNMoments();
    _filterFunction   = Vector( _settings->GetNTotal(), 1.0 );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index =
                    unsigned( ( i - i % unsigned( std::pow( nMoments + 1, l ) ) ) / unsigned( std::pow( nMoments + 1, l ) ) ) % ( nMoments + 1 );
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

        double thetaTilde1 = pow( rho[k], 2 ) * pow( ETilde, 2 ) - 2 * rho[k] * rhoTilde * E[k] * ETilde + 2 * rho[k] * E[k] * pow( m1Tilde, 2 ) +
                             2 * rho[k] * E[k] * pow( m2Tilde, 2 ) - 2 * rho[k] * ETilde * m1[k] * m1Tilde - 2 * rho[k] * ETilde * m2[k] * m2Tilde;

        double thetaTilde2 = pow( rhoTilde, 2 ) * pow( E[k], 2 ) - 2 * rhoTilde * E[k] * m1[k] * m1Tilde - 2 * rhoTilde * E[k] * m2[k] * m2Tilde +
                             2 * rhoTilde * ETilde * pow( m1[k], 2 ) + 2 * rhoTilde * ETilde * pow( m2[k], 2 ) -
                             pow( m1[k] * m2Tilde - m2[k] * m1Tilde, 2 );

        double thetaTilde = thetaTilde1 + thetaTilde2;

        double denominator = pow( m1[k] - m1Tilde, 2 ) + pow( m2[k] - m2Tilde, 2 ) - 2 * rho[k] * E[k] + 2 * rho[k] * ETilde + 2 * rhoTilde * E[k] -
                             2 * rhoTilde * ETilde;

        double part1Nominator = rho[k] * ETilde - 2 * rho[k] * E[k] + rhoTilde * E[k] - m1[k] * m1Tilde - m2[k] * m2Tilde;

        double t2a = ( part1Nominator + sqrt( thetaTilde ) ) / denominator;
        double t2b = ( part1Nominator - sqrt( thetaTilde ) ) / denominator;

        if( t2a > 1.0 || t2a < 0.0 || !std::isfinite( t2a ) ) t2a = 0.0;
        if( t2b > 1.0 || t2b < 0.0 || !std::isfinite( t2b ) ) t2b = 0.0;

        double t2MaxTmp = MathTools::max( t2a, t2b );
        if( t2MaxTmp > t2Max ) t2Max = t2MaxTmp;
    }

    double theta = MathTools::max( t2Max, t1Max );
    std::cout << t2Max << std::endl;
    std::cout << "theta = " << theta << std::endl;

    lambda = theta * u2 + ( 1 - theta ) * u;
}

void HyperbolicityLimiter2D::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) { this->SolveClosure( lambda, u, refLevel ); }
