#include "filters/fokkerplanckfilter.h"

FokkerPlanckFilter::FokkerPlanckFilter( Settings* settings ) : Filter( settings ) {
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[FokkerPlanckFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

void FokkerPlanckFilter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, m ) ) ) / unsigned( std::pow( maxDegree + 1, m ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= std::exp( -_lambda * index * ( index + 1 ) );
        }
    }
    _lambda = 0.0;    // reset lambda
}

void FokkerPlanckFilter::FilterMoments( Tensor& u ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                u( s, l, i ) = _filterFunction[i] * u( s, l, i );
            }
        }
    }
}

void FokkerPlanckFilter::FilterMoments( Tensor& v, const Tensor& u ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                v( s, l, i ) = _filterFunction[i] * u( s, l, i );
            }
        }
    }
}

void FokkerPlanckFilter::FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < v.columns(); ++i ) {
            v( s, i ) = _filterFunction[i] * u( s, l, i );
        }
    }
}

FokkerPlanckFilter::~FokkerPlanckFilter() {}
