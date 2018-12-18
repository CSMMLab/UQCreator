#include "lassofilter.h"

LassoFilter::LassoFilter( Settings* settings ) : Closure( settings ) {
    _alpha            = 1.0;    // unsigned n;
    unsigned nMoments = _settings->GetNMoments();
    _filterParam      = Vector( _settings->GetNTotal(), 1.0 );
    _l1Norms          = Vector( _settings->GetNTotal(), 0.0 );

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned k = 0; k < _settings->GetNQTotal(); ++k ) {
            _l1Norms[i] += std::fabs( _phiTildeWf( k, i ) );
        }
    }

    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index = unsigned( ( i - i % unsigned( std::pow( nMoments, l ) ) ) / unsigned( std::pow( nMoments, l ) ) ) % nMoments;
            _filterParam[i] *= index * ( index + 1 );
        }
    }
}

LassoFilter::~LassoFilter() {}

void LassoFilter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void LassoFilter::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix LassoFilter::U( const Matrix& Lambda ) { return Lambda; }

void LassoFilter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void LassoFilter::SolveClosure( Matrix& lambda, const Matrix& u ) {
    double scL1, filterStrength;
    unsigned nMax = _settings->GetNTotal() - 1;

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            filterStrength = std::fabs( u( s, nMax ) ) / ( _filterParam[nMax] * _l1Norms[nMax] );
            scL1           = 1.0 - filterStrength * _filterParam[i] * _l1Norms[i] / std::fabs( u( s, i ) );
            if( scL1 < 0 || std::fabs( u( s, i ) ) < 1e-7 ) scL1 = 0.0;
            lambda( s, i ) = scL1 * u( s, i );
        }
    }
}
