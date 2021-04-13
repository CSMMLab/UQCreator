#include "filters/houlifilter.h"

HouLiFilter::HouLiFilter( Settings* settings ) : Filter( settings ) {
    _gamma = 36;
    _eps   = 1.0 / _lambda;

    // reset lambda to 1
    _lambda = 1.0;

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[HouLiFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

void HouLiFilter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= exp( -FilterFunction( double( index ) / double( maxDegree + 1 ) ) / _eps );
        }
    }
}

HouLiFilter::~HouLiFilter() {}

double HouLiFilter::FilterFunction( double eta ) const {
    if( eta <= 2.0 / 3.0 ) {
        return 0.0;
    }
    else {
        return pow( eta, _gamma );
    }
}
