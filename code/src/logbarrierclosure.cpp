#include "logbarrierclosure.h"

LogBarrierClosure::LogBarrierClosure( Settings* settings ) : Closure( settings ) {
    _alpha    = 1.0;
    double du = 0.5;
    // TODO: IC1, IC3 Test cases should not be hard coded
    if( _settings->GetNDimXi() == 2 ) {
        _uMinus = 1.0 - du;
        _uPlus  = 12.0 + du + 0.2;
    }
    else {
        _uMinus = 3.0 - du;
        _uPlus  = 12.0 + du;
    }
}

LogBarrierClosure::~LogBarrierClosure() {}

void LogBarrierClosure::U( Vector& out, const Vector& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        out[l] = 0.5 * ( _uMinus + _uPlus +
                         pow( _uMinus - _uPlus, 2 ) * Lambda[l] / ( sqrt( pow( _uMinus - _uPlus, 2 ) * pow( Lambda[l], 2 ) + 4.0 ) + 2.0 ) );
    }
}

void LogBarrierClosure::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
            out( l, k ) =
                0.5 * ( _uMinus + _uPlus +
                        pow( _uMinus - _uPlus, 2 ) * Lambda( l, k ) / ( sqrt( pow( _uMinus - _uPlus, 2 ) * pow( Lambda( l, k ), 2 ) + 4.0 ) + 2.0 ) );
        }
    }
}

Matrix LogBarrierClosure::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
            y( l, k ) =
                0.5 * ( _uMinus + _uPlus +
                        pow( _uMinus - _uPlus, 2 ) * Lambda( l, k ) / ( sqrt( pow( _uMinus - _uPlus, 2 ) * pow( Lambda( l, k ), 2 ) + 4.0 ) + 2.0 ) );
        }
    }
    return y;
}

double LogBarrierClosure::U( const double Lambda ) {
    return 0.5 * ( _uMinus + _uPlus + pow( _uMinus - _uPlus, 2 ) * Lambda / ( sqrt( pow( _uMinus - _uPlus, 2 ) * pow( Lambda, 2 ) + 4.0 ) + 2.0 ) );
}

void LogBarrierClosure::DU( Matrix& y, const Vector& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned m = 0; m < _nStates; ++m ) {
            y( l, m ) = 1.0 / ( 1.0 / ( pow( U( Lambda[l] ) - _uMinus, 2 ) ) + 1. / ( pow( _uPlus - U( Lambda[l] ), 2 ) ) );
        }
    }
}
