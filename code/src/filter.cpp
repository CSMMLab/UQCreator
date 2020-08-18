#include "filter.h"
#include "erfcfilter.h"
#include "erfclogfilter.h"
#include "exponentialfilter.h"
#include "fokkerplanckfilter.h"
#include "houlifilter.h"
#include "l2filter.h"
#include "lassofilter.h"
#include "splinefilter.h"

Filter::Filter( Settings* settings ) : _settings( settings ), _lambda( _settings->GetFilterStrength() ) {
    _log            = spdlog::get( "event" );
    _filterFunction = Vector( _settings->GetNTotal(), 1.0 );

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Filter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Filter::~Filter() {}

void Filter::SetupFilter() {
    // turn off filter
    _lambda = 0.0;
}

void Filter::FilterMoments( Tensor& u ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                u( s, l, i ) = pow( _filterFunction[i], _settings->GetDT() * _lambda ) * u( s, l, i );
            }
        }
    }
}

void Filter::FilterMoments( Tensor& v, const Tensor& u ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                v( s, l, i ) = pow( _filterFunction[i], _settings->GetDT() * _lambda ) * u( s, l, i );
            }
        }
    }
}

void Filter::FilterMoments( Matrix& v, const Tensor& u, unsigned l ) const {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < v.columns(); ++i ) {
            v( s, i ) = pow( _filterFunction[i], _settings->GetDT() * _lambda ) * u( s, l, i );
        }
    }
}

Filter* Filter::Create( Settings* settings ) {
    auto log        = spdlog::get( "event" );
    auto filterType = settings->GetFilterType();
    Filter* tmp;
    if( filterType == FilterType::F_L2FILTER ) {
        tmp = new L2Filter( settings );
    }
    else if( filterType == FilterType::F_LASSOFILTER ) {
        tmp = new LassoFilter( settings );
    }
    else if( filterType == FilterType::F_EXPFILTER ) {
        tmp = new ExponentialFilter( settings );
    }
    else if( filterType == FilterType::F_SPLINEFILTER ) {
        tmp = new SplineFilter( settings );
    }
    else if( filterType == FilterType::F_HOULIFILTER ) {
        tmp = new HouLiFilter( settings );
    }
    else if( filterType == FilterType::F_FOKKERPLANCKFILTER ) {
        tmp = new FokkerPlanckFilter( settings );
    }
    else if( filterType == FilterType::F_ERFCFILTER ) {
        tmp = new ErfcFilter( settings );
    }
    else if( filterType == FilterType::F_ERFCLOGFILTER ) {
        tmp = new ErfcLogFilter( settings );
    }
    else if( filterType == FilterType::F_NOFILTER ) {
        tmp = new Filter( settings );
    }
    else {
        log->error( "[filter]: Invalid filter type" );
        exit( EXIT_FAILURE );
    }
    tmp->SetupFilter();
    return tmp;
}
