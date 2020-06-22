#include "kineticclosure.h"

KineticClosure::KineticClosure( Settings* settings ) : Closure( settings ), _nMoments( 2 ) { _alpha = 1.0; }

KineticClosure::~KineticClosure() {}

void KineticClosure::U( Vector& out, const Vector& Lambda ) { out[0] = exp( Lambda[0] ); }

void KineticClosure::U( Vector& out, const Vector& Lambda, bool dummy ) { out[0] = exp( Lambda[0] ); }

void KineticClosure::U( Tensor& out, const Tensor& Lambda ) {
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            out( 0, l, k ) = exp( Lambda( 0, l, k ) );
        }
    }
}

Tensor KineticClosure::U( const Tensor& Lambda ) {
    Tensor y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            y( 0, l, k ) = exp( Lambda( 0, l, k ) );
        }
    }

    return y;
}

void KineticClosure::DU( Matrix& y, const Vector& Lambda ) { y( 0, 0 ) = exp( Lambda[0] ); }
