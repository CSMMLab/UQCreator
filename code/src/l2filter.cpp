#include "l2filter.h"

L2Filter::L2Filter( Settings* settings ) : Closure( settings ), _lambda( 0.00001 ) { _alpha = 1.0; }

L2Filter::~L2Filter() {}

void L2Filter::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void L2Filter::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix L2Filter::U( const Matrix& Lambda ) { return Lambda; }

void L2Filter::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void L2Filter::SolveClosure( Matrix& lambda, const Matrix& u ) {
    // unsigned n;
    unsigned nMoments = _settings->GetNMoments();
    double tmp;
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            tmp = 1.0;
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index = unsigned( ( i - i % unsigned( std::pow( nMoments, l ) ) ) / unsigned( std::pow( nMoments, l ) ) ) % nMoments;
                tmp *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) );
            }
            lambda( s, i ) = tmp * u( s, i );
        }
    }
}
