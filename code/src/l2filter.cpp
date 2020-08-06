#include "l2filter.h"

L2Filter::L2Filter( Settings* settings ) : Filter( settings ) {}

void L2Filter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
        for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
            double eta = _lambda * pow( 1, 2 ) * pow( 1 + 1, 2 );    // punishes variance dampening
            // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            unsigned index =
                unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
            if( i == 0 )
                _filterFunction[i] *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) );
            else
                _filterFunction[i] *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) - eta );
        }
    }
    _lambda = 0.0;    // reset lambda
}

L2Filter::~L2Filter() {}
