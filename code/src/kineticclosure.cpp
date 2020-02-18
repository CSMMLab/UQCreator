#include "kineticclosure.h"

KineticClosure::KineticClosure( Settings* settings ) : Closure( settings ), _nMoments( 2 ) { _alpha = 1.0; }

KineticClosure::~KineticClosure() {}

void KineticClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = exp( Lambda[0] );
}

void KineticClosure::U( Vector& out, const Vector& Lambda, bool dummy ) { out[0] = exp( Lambda[0] ); }

void KineticClosure::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        out( 0, k ) = exp( Lambda( 0, k ) );
    }
}

Matrix KineticClosure::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        y( 0, k ) = exp( Lambda( 0, k ) );
    }

    return y;
}

void KineticClosure::DU( Matrix& y, const Vector& Lambda ) { y( 0, 0 ) = exp( Lambda[0] ); }
