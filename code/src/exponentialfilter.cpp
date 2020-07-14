#include "exponentialfilter.h"
#include <limits>

ExponentialFilter::ExponentialFilter( Settings* settings ) : Closure( settings ), _lambda( _settings->GetFilterStrength() ) {
    _alpha             = 1.0;    // unsigned n;
    double epsilonM    = std::numeric_limits<double>::denorm_min();
    _c                 = log( epsilonM );
    unsigned maxDegree = _settings->GetMaxDegree();
    _filterFunction    = Vector( _settings->GetNTotal(), 1.0 );

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ExponentialFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= pow( FilterFunction( double( index ) / double( maxDegree + 1 ) ), _lambda );
        }
    }
}

ExponentialFilter::~ExponentialFilter() {}

double ExponentialFilter::FilterFunction( double eta ) const { return exp( _c * pow( eta, _filterOrder ) ); }

void ExponentialFilter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void ExponentialFilter::U( Tensor& out, const Tensor& Lambda ) { out = Lambda; }

Tensor ExponentialFilter::U( const Tensor& Tensor ) { return Tensor; }

void ExponentialFilter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void ExponentialFilter::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                lambda( s, l, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, l, i );
            }
        }
    }
}

void ExponentialFilter::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                lambda( s, l, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, l, i );
            }
        }
    }
}
