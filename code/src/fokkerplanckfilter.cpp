#include "fokkerplanckfilter.h"

FokkerPlanckFilter::FokkerPlanckFilter( Settings* settings ) : Filter( settings ) {}

void FokkerPlanckFilter::SetupFilter() {
    unsigned maxDegree = _settings->GetMaxDegree();
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                unsigned index =
                    unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
                _filterFunction[i] *= std::exp( -_lambda * index * ( index + 1 ) );
            }
        }
    }
    _lambda = 0.0;    // reset lambda
}

FokkerPlanckFilter::~FokkerPlanckFilter() {}
