#include "logsin.h"

LogSin::LogSin( Settings* settings ) : Closure( settings ) {
    _alpha    = 1.0;
    double du = 0.1;
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

LogSin::~LogSin() {}

void LogSin::U( Vector& out, const Vector& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        out[l] = ( ( _uPlus - _uMinus ) / PI ) * ( PI / 2 + std::atan( Lambda[l] ) ) + _uMinus;
    }
}

void LogSin::U( Tensor& out, const Tensor& Lambda ) {
    for( unsigned n = 0; n < _nMultiElements; ++n ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
                out( l, n, k ) = ( ( _uPlus - _uMinus ) / PI ) * ( PI / 2 + std::atan( Lambda( l, n, k ) ) ) + _uMinus;
            }
        }
    }
}

Tensor LogSin::U( const Tensor& Lambda ) {
    Tensor y( _nStates, _nMultiElements, Lambda.columns(), 0.0 );
    for( unsigned n = 0; n < _nMultiElements; ++n ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned int k = 0; k < Lambda.columns(); ++k ) {
                y( n, l, k ) = ( ( _uPlus - _uMinus ) / PI ) * ( PI / 2 + std::atan( Lambda( n, l, k ) ) ) + _uMinus;
            }
        }
    }
    return y;
}

void LogSin::DU( Matrix& y, const Vector& Lambda ) {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned m = 0; m < _nStates; ++m ) {
            y( l, m ) = ( ( _uPlus - _uMinus ) / PI ) / ( std::pow( Lambda[l], 2 ) + 1.0 );
        }
    }
}
