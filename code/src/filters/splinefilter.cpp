#include "filters/splinefilter.h"

SplineFilter::SplineFilter( Settings* settings ) : Filter( settings ) {
    unsigned maxDegree = _settings->GetMaxDegree();
    _eta               = pow( 1.0 / double( maxDegree + 1 ), 4 );    // set to zero to turn off variance correction

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[SplineFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

void SplineFilter::SetupFilter() {
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

SplineFilter::~SplineFilter() {}

double SplineFilter::FilterFunction( double eta ) const {
    if( eta > 0.1 / double( _settings->GetMaxDegree() + 1 ) )
        return 1.0 / ( pow( eta, 4 ) + 1.0 - _eta );
    else
        return 1.0 / ( pow( eta, 4 ) + 1.0 );
}
