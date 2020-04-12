#include "houlifilter.h"

HouLiFilter::HouLiFilter( Settings* settings ) : Closure( settings ), _lambda( _settings->GetFilterStrength() ) {
    _alpha            = 1.0;    // unsigned n;
    double epsilonM   = std::numeric_limits<double>::denorm_min();
    _gamma            = 36;
    double _eps       = 1.0 / _lambda;
    unsigned nMoments = _settings->GetNMoments();
    _filterFunction   = Vector( _settings->GetNTotal(), 1.0 );

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "moment_system" );
        _filterOrder = problem->get_as<double>( "filterOrder" ).value_or( 1 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[HouLiFilter] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( nMoments + 1, l ) ) ) / unsigned( std::pow( nMoments + 1, l ) ) ) % ( nMoments + 1 );
            _filterFunction[i] *= exp( -FilterFunction( double( index ) / double( nMoments + 1 ) ) / _eps );
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

void HouLiFilter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void HouLiFilter::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix HouLiFilter::U( const Matrix& Lambda ) { return Lambda; }

void HouLiFilter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void HouLiFilter::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, i );
        }
    }
}

void HouLiFilter::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = pow( _filterFunction[i], _settings->GetDT() ) * u( s, i );
        }
    }
}
