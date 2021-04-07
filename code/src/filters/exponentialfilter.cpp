#include "filters/exponentialfilter.h"
#include <limits>

ExponentialFilter::ExponentialFilter( Settings* settings ) : Filter( settings ) {
    double epsilonM = std::numeric_limits<double>::denorm_min();
    _c              = log( epsilonM );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ExponentialFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ExponentialFilter::~ExponentialFilter() {}

void ExponentialFilter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= FilterFunction( double( index ) / double( maxDegree + 1 ) );
        }
    }
}

double ExponentialFilter::FilterFunction( double eta ) const { return exp( _c * pow( eta, _filterOrder ) ); }
