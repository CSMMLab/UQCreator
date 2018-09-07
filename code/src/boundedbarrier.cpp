#include "boundedbarrier.h"

BoundedBarrier::BoundedBarrier( Problem* problem ) : Closure( problem ) {}

BoundedBarrier::~BoundedBarrier() {}

void BoundedBarrier::U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda ) {
    double ePos, eNeg;
    for( int l = 0; l < _nStates; ++l ) {
        ePos = exp( Lambda[l] );
        eNeg = 1.0 / ePos;
        if( Lambda[l] > 0.0 ) {
            out[l] = _uPlus / ( eNeg + 1.0 ) + _uMinus * eNeg / ( 1.0 + eNeg );
        }
        else {
            out[l] = _uMinus / ( ePos + 1.0 ) + _uPlus * ePos / ( 1.0 + ePos );
        }
    }
}

blaze::DynamicMatrix<double> BoundedBarrier::U( const blaze::DynamicMatrix<double>& Lambda ) {
    double ePos, eNeg;
    blaze::DynamicMatrix<double> y( _nStates, Lambda.columns(), 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
            ePos = exp( Lambda( l, k ) );
            eNeg = 1.0 / ePos;
            if( Lambda( l, k ) > 0 ) {
                y( l, k ) = _uPlus / ( eNeg + 1.0 ) + _uMinus * eNeg / ( 1.0 + eNeg );
            }
            else {
                y( l, k ) = _uMinus / ( ePos + 1.0 ) + _uPlus * ePos / ( 1.0 + ePos );
            }
        }
    }
    return y;
}

void BoundedBarrier::DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda ) {
    for( int l = 0; l < _nStates; ++l ) {
        for( int m = 0; m < _nStates; ++m ) {
            double ePos = exp( Lambda[l] );
            double eNeg = 1 / ePos;
            if( Lambda[l] > 0 ) {
                y( l, m ) = eNeg * ( _uPlus - _uMinus ) / ( exp( -2.0 * Lambda[l] ) + 2.0 * eNeg + 1.0 );
            }
            else {
                y( l, m ) = ePos * ( _uPlus - _uMinus ) / ( 1.0 + 2.0 * ePos + exp( 2.0 * Lambda[l] ) );
            }
        }
    }
}
