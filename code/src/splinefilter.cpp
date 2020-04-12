#include "splinefilter.h"

SplineFilter::SplineFilter( Settings* settings ) : Closure( settings ), _lambda( _settings->GetFilterStrength() ) {
    _alpha            = 1.0;    // unsigned n;
    unsigned nMoments = _settings->GetNMoments();
    _filterFunction   = Vector( _settings->GetNTotal(), 1.0 );
    _eta              = pow( 1.0 / double( nMoments + 1 ), 4 );    // set to zero to turn off variance correction

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
                unsigned( ( i - i % unsigned( std::pow( nMoments + 1, l ) ) ) / unsigned( std::pow( nMoments + 1, l ) ) ) % ( nMoments + 1 );
            _filterFunction[i] *= pow( FilterFunction( double( index ) / double( nMoments + 1 ) ), _lambda );
        }
        std::cout << "g = " << _filterFunction[i] << std::endl;
    }
    std::cout << "-------------" << std::endl;
    // exit( EXIT_FAILURE );
}

SplineFilter::~SplineFilter() {}

double SplineFilter::FilterFunction( double eta ) const {
    if( eta > 0.1 / double( _settings->GetNMoments() + 1 ) )
        return 1.0 / ( pow( eta, 4 ) + 1.0 - _eta );
    else
        return 1.0 / ( pow( eta, 4 ) + 1.0 );
}

void SplineFilter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void SplineFilter::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix SplineFilter::U( const Matrix& Lambda ) { return Lambda; }

void SplineFilter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void SplineFilter::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, i );
        }
    }
}

void SplineFilter::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, i );
        }
    }
}
