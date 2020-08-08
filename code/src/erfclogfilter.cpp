#include "erfclogfilter.h"

ErfcLogFilter::ErfcLogFilter( Settings* settings ) : Filter( settings ) {
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ErfcLogFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ErfcLogFilter::~ErfcLogFilter() {}

void ErfcLogFilter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            _filterFunction[i] *= FilterFunction( double( index ) / double( maxDegree + 1 ) );
        }
        std::cout << _filterFunction[i] << std::endl;
    }
}

double ErfcLogFilter::FilterFunction( double eta ) const {
    double thetaBar = std::fabs( eta ) - 0.5;
    return 0.5 * erfc( 2.0 * sqrt( _filterOrder ) * thetaBar * sqrt( -log( 1.0 - 4.0 * pow( thetaBar, 2 ) ) / ( 4.0 * pow( thetaBar, 2 ) ) ) );
}
