#include "l2filter.h"

L2Filter::L2Filter( Settings* settings ) : Closure( settings ), _lambda( 0.0000015 ) {
    _alpha            = 1.0;    // unsigned n;
    unsigned nMoments = _settings->GetNMoments();
    _filterFunction   = Vector( _settings->GetNTotal(), 1.0 );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index = unsigned( ( i - i % unsigned( std::pow( nMoments, l ) ) ) / unsigned( std::pow( nMoments, l ) ) ) % nMoments;
                _filterFunction[i] *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) );
            }
        }
    }
}

L2Filter::~L2Filter() {}

void L2Filter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void L2Filter::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix L2Filter::U( const Matrix& Lambda ) { return Lambda; }

void L2Filter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void L2Filter::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = _filterFunction[i] * u( s, i );
        }
    }
}

void L2Filter::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            lambda( s, i ) = _filterFunction[i] * u( s, i );
        }
    }
}
