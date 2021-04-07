#include "closures/boundedbarrier.h"

BoundedBarrier::BoundedBarrier( Settings* settings ) : Closure( settings ) {
    _alpha    = 1.0;
    double du = 0.0001;
    // TODO: IC1, IC3 Test cases should not be hard coded
    if( _settings->GetNDimXi() == 2 ) {
        _uMinus = 1.0 - du;
        _uPlus  = 12.0 + du + 0.2;
    }
    else {
        _uMinus = 1.0 - du;
        _uPlus  = 12.0 + du;
    }
}

BoundedBarrier::~BoundedBarrier() {}

void BoundedBarrier::U( Vector& out, const Vector& Lambda ) {
    double ePos, eNeg;
    for( unsigned l = 0; l < _nStates; ++l ) {
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

void BoundedBarrier::U( Tensor& out, const Tensor& Lambda ) {
    double ePos, eNeg;
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
                ePos = exp( Lambda( s, l, k ) );
                eNeg = 1.0 / ePos;
                if( Lambda( s, l, k ) > 0 ) {
                    out( s, l, k ) = _uPlus / ( eNeg + 1.0 ) + _uMinus * eNeg / ( 1.0 + eNeg );
                }
                else {
                    out( s, l, k ) = _uMinus / ( ePos + 1.0 ) + _uPlus * ePos / ( 1.0 + ePos );
                }
            }
        }
    }
}

Tensor BoundedBarrier::U( const Tensor& Lambda ) {
    double ePos, eNeg;
    Tensor y( _nStates, _nMultiElements, Lambda.columns(), 0.0 );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
                ePos = exp( Lambda( s, l, k ) );
                eNeg = 1.0 / ePos;
                if( Lambda( s, l, k ) > 0 ) {
                    y( s, l, k ) = _uPlus / ( eNeg + 1.0 ) + _uMinus * eNeg / ( 1.0 + eNeg );
                }
                else {
                    y( s, l, k ) = _uMinus / ( ePos + 1.0 ) + _uPlus * ePos / ( 1.0 + ePos );
                }
            }
        }
    }
    return y;
}

void BoundedBarrier::DU( Matrix& y, const Vector& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned m = 0; m < _nStates; ++m ) {
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
