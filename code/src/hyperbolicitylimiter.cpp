#include "hyperbolicitylimiter.h"
#include "mathtools.h"

HyperbolicityLimiter::HyperbolicityLimiter( Settings* settings )
    : Closure( settings ), _gamma( _settings->GetGamma() ), _lambda( _settings->GetFilterStrength() ) {
    _alpha          = 1.0;
    double epsilonM = std::numeric_limits<double>::denorm_min();
    _c              = log( epsilonM );
    if( _lambda < 0 ) _lambda = 0.0;
    unsigned maxDegree = _settings->GetMaxDegree();
    _filterFunction    = Vector( _settings->GetNTotal(), 1.0 );

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[SplineFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }

    std::cout << "-------------" << std::endl;
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= pow( FilterFunction( double( index ) / double( maxDegree + 1 ) ), _lambda );
        }
        std::cout << "g = " << _filterFunction[i] << std::endl;
    }
    std::cout << "-------------" << std::endl;
    // exit( EXIT_FAILURE );
}

HyperbolicityLimiter::~HyperbolicityLimiter() {}

double HyperbolicityLimiter::FilterFunction( double eta ) const { return exp( _c * pow( eta, _filterOrder ) ); }

void HyperbolicityLimiter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void HyperbolicityLimiter::U( Tensor& out, const Tensor& Lambda ) { out = Lambda; }

Tensor HyperbolicityLimiter::U( const Tensor& Lambda ) { return Lambda; }

void HyperbolicityLimiter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void HyperbolicityLimiter::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    Vector rho( _nQTotalForRef[refLevel], 0.0 );
    Vector m( _nQTotalForRef[refLevel], 0.0 );
    Vector E( _nQTotalForRef[refLevel], 0.0 );
    Vector p( _nQTotalForRef[refLevel], 0.0 );
    Vector t1( _nQTotalForRef[refLevel], 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];

    for( unsigned l = 0; l < _nMultiElements; ++l ) {

        double t1Max = -1000;
        double t2Max = -1000;

        // save zero order moments (realizable solution)
        double rhoTilde = u( 0, l, 0 );
        double mTilde   = u( 1, l, 0 );
        double ETilde   = u( 2, l, 0 );

        // save zero order moment as full moment vector for output
        Matrix u2( _nStates, nTotal, 0.0 );
        u2( 0, 0 ) = u( 0, l, 0 );
        u2( 1, 0 ) = u( 1, l, 0 );
        u2( 2, 0 ) = u( 2, l, 0 );

        double pTilde = ( _gamma - 1.0 ) * ( ETilde - 0.5 * pow( mTilde, 2 ) / rhoTilde );

        // Question: Is e a scalar in our case?
        double e = MathTools::min( 1e-10, MathTools::min( rhoTilde, pTilde ) );

        // determine t1Max
        for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
            Vector uKinetic = EvaluateLambda( u, l, k, nTotal );

            rho[k] = uKinetic[0];
            m[k]   = uKinetic[1];
            E[k]   = uKinetic[2];

            t1[k] = ( e - rho[k] ) / ( rhoTilde - rho[k] );
            // t1[k] = rho[k] / ( rho[k] - rhoTilde );
            if( t1[k] > 1.0 || t1[k] < 0.0 || !std::isfinite( t1[k] ) ) t1[k] = 0.0;

            if( t1[k] > t1Max ) t1Max = t1[k];
        }

        // determine t2Max
        for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
            // rho[k] = t1Max * rhoTilde + ( 1.0 - t1Max ) * rho[k];

            // p[k] = ( _gamma - 1.0 ) * ( E[k] - 0.5 * pow( m[k], 2 ) / rho[k] );

            // a* t ^ 2 + b* t + c = 0
            double a = ( E[k] - ETilde ) * ( rho[k] - rhoTilde ) - 0.5 * pow( m[k] - mTilde, 2 );
            double b = m[k] * ( m[k] - mTilde ) - rho[k] * ( E[k] - ETilde ) - E[k] * ( rho[k] - rhoTilde );
            // double b = m[k] * ( m[k] - mTilde ) + rho[k] * ETilde + E[k] * rhoTilde - 2.0 * E[k] * rho[k];    // my calculation of b is the same
            double c = E[k] * rho[k] - 0.5 * pow( m[k], 2 ) - e;

            // Citardauq Formula as stable p-q formula
            double q = -0.5 * ( b + MathTools::sign( b ) * sqrt( pow( b, 2 ) - 4.0 * a * c ) );
            if( pow( b, 2 ) - 4.0 * a * c < 0 ) std::cerr << "[RealizabilityLimiter]: Imaginary unit detected!" << std::endl;
            double t2a = q / a;
            double t2b = c / q;
            /*
                        // test root
                        double rhoA = t2a * rhoTilde + ( 1.0 - t2a ) * rho[k];
                        double Ea   = t2a * ETilde + ( 1.0 - t2a ) * E[k];
                        double ma   = t2a * mTilde + ( 1.0 - t2a ) * m[k];
                        double root = rhoA * Ea - 0.5 * pow( ma, 2 );
                        // std::cout << "theta = " << t2a << ":   root A " << root << std::endl;

                        double rhoB = t2b * rhoTilde + ( 1.0 - t2b ) * rho[k];
                        double Eb   = t2b * ETilde + ( 1.0 - t2b ) * E[k];
                        double mb   = t2b * mTilde + ( 1.0 - t2b ) * m[k];
                        root        = rhoB * Eb - 0.5 * pow( mb, 2 );
                        // std::cout << "theta = " << t2b << "root B " << root << std::endl;
            */
            if( t2a > 1.0 || t2a < 0.0 || !std::isfinite( t2a ) ) t2a = 0.0;
            if( t2b > 1.0 || t2b < 0.0 || !std::isfinite( t2b ) ) t2b = 0.0;

            double t2MaxTmp = MathTools::max( t2a, t2b );
            if( t2MaxTmp > t2Max ) t2Max = t2MaxTmp;
        }
        double theta = t2Max;
        if( t1Max > t2Max ) theta = t1Max;
        theta = MathTools::max( theta, 0.0 );
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                lambda( s, l, i ) = theta * u2( s, i ) + ( 1.0 - theta ) * u( s, l, i );
            }
        }
    }
}

void HyperbolicityLimiter::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    Vector rho( _nQTotalForRef[refLevel], 0.0 );
    Vector m( _nQTotalForRef[refLevel], 0.0 );
    Vector E( _nQTotalForRef[refLevel], 0.0 );
    Vector p( _nQTotalForRef[refLevel], 0.0 );
    Vector t1( _nQTotalForRef[refLevel], 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    Tensor uF       = u;

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                uF( s, l, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, l, i );
            }
        }
    }

    for( unsigned l = 0; l < _nMultiElements; ++l ) {

        double t1Max = -1000;
        double t2Max = -1000;

        // save zero order moments (realizable solution)
        double rhoTilde = uF( 0, l, 0 );
        double mTilde   = uF( 1, l, 0 );
        double ETilde   = uF( 2, l, 0 );

        // save zero order moment as full moment vector for output
        Matrix u2( _nStates, nTotal, 0.0 );
        u2( 0, 0 ) = uF( 0, l, 0 );
        u2( 1, 0 ) = uF( 1, l, 0 );
        u2( 2, 0 ) = uF( 2, l, 0 );

        double pTilde = ( _gamma - 1.0 ) * ( ETilde - 0.5 * pow( mTilde, 2 ) / rhoTilde );

        // Question: Is e a scalar in our case?
        double e = MathTools::min( 1e-10, MathTools::min( rhoTilde, pTilde ) );

        // determine t1Max
        for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
            Vector uKinetic = EvaluateLambda( uF, l, k, nTotal );

            rho[k] = uKinetic[0];
            m[k]   = uKinetic[1];
            E[k]   = uKinetic[2];

            t1[k] = ( e - rho[k] ) / ( rhoTilde - rho[k] );
            // t1[k] = rho[k] / ( rho[k] - rhoTilde );
            if( t1[k] > 1.0 || t1[k] < 0.0 || !std::isfinite( t1[k] ) ) t1[k] = 0.0;

            if( t1[k] > t1Max ) t1Max = t1[k];
        }

        // determine t2Max
        for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
            // rho[k] = t1Max * rhoTilde + ( 1.0 - t1Max ) * rho[k];

            // p[k] = ( _gamma - 1.0 ) * ( E[k] - 0.5 * pow( m[k], 2 ) / rho[k] );

            // a* t ^ 2 + b* t + c = 0
            double a = ( E[k] - ETilde ) * ( rho[k] - rhoTilde ) - 0.5 * pow( m[k] - mTilde, 2 );
            double b = m[k] * ( m[k] - mTilde ) - rho[k] * ( E[k] - ETilde ) - E[k] * ( rho[k] - rhoTilde );
            // double b = m[k] * ( m[k] - mTilde ) + rho[k] * ETilde + E[k] * rhoTilde - 2.0 * E[k] * rho[k];    // my calculation of b is the same
            double c = E[k] * rho[k] - 0.5 * pow( m[k], 2 ) - e;

            // Citardauq Formula as stable p-q formula
            double q = -0.5 * ( b + MathTools::sign( b ) * sqrt( pow( b, 2 ) - 4.0 * a * c ) );
            if( pow( b, 2 ) - 4.0 * a * c < 0 ) std::cerr << "[RealizabilityLimiter]: Imaginary unit detected!" << std::endl;
            double t2a = q / a;
            double t2b = c / q;
            /*
                        // test root
                        double rhoA = t2a * rhoTilde + ( 1.0 - t2a ) * rho[k];
                        double Ea   = t2a * ETilde + ( 1.0 - t2a ) * E[k];
                        double ma   = t2a * mTilde + ( 1.0 - t2a ) * m[k];
                        double root = rhoA * Ea - 0.5 * pow( ma, 2 );
                        // std::cout << "theta = " << t2a << ":   root A " << root << std::endl;

                        double rhoB = t2b * rhoTilde + ( 1.0 - t2b ) * rho[k];
                        double Eb   = t2b * ETilde + ( 1.0 - t2b ) * E[k];
                        double mb   = t2b * mTilde + ( 1.0 - t2b ) * m[k];
                        root        = rhoB * Eb - 0.5 * pow( mb, 2 );
                        // std::cout << "theta = " << t2b << "root B " << root << std::endl;
                        */

            if( t2a > 1.0 || t2a < 0.0 || !std::isfinite( t2a ) ) t2a = 0.0;
            if( t2b > 1.0 || t2b < 0.0 || !std::isfinite( t2b ) ) t2b = 0.0;

            double t2MaxTmp = MathTools::max( t2a, t2b );
            if( t2MaxTmp > t2Max ) t2Max = t2MaxTmp;
        }
        double theta = t2Max;
        if( t1Max > t2Max ) theta = t1Max;
        theta = MathTools::max( theta, 0.0 );
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                lambda( s, l, i ) = theta * u2( s, i ) + ( 1.0 - theta ) * uF( s, l, i );
            }
        }
    }

    /*
        // test positivity for limited moments
        for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
            Vector uKinetic = EvaluateLambda( lambda,l, k, nTotal );

            rho[k] = uKinetic[0];
            m[k]   = uKinetic[1];
            E[k]   = uKinetic[2];
            p[k]   = ( _gamma - 1.0 ) * ( E[k] - 0.5 * pow( m[k], 2 ) / rho[k] );

            if( rho[k] < 0.0 || E[k] < 0.0 || p[k] < 0.0  ) {
                std::cerr << "ERROR: Non-realizable reconstruction detected" << std::endl;
                std::cout << "theta = " << theta << std::endl;
                std::cout << "Moment vectors " << uF << std::endl;
                std::cout << "states are " << rho[k] << " " << E[k] << " " << p[k] << std::endl;
            }
            if( theta > 1e-5 ) {
                std::cout << "states are " << rho[k] << " " << E[k] << " " << p[k] << " - theta = " << theta << std::endl;
            }
        }*/
}
