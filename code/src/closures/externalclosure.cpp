#include "closures/externalclosure.h"

ExternalClosure::ExternalClosure( void ( *U )( double*, double* ),
                                  void ( *DU )( double*, double* ),
                                  void ( *SolveClosure )( double*, double*, unsigned ),
                                  Settings* settings )
    : Closure( settings ), _U( U ), _DU( DU ), _SolveClosure( SolveClosure ) {}

ExternalClosure::~ExternalClosure() {}

void ExternalClosure::U( Vector& /*out*/, const Vector& /*Lambda*/ ) { exit( EXIT_FAILURE ); }

void ExternalClosure::U( Tensor& out, const Tensor& Lambda ) { _U( out.GetPointer(), Lambda.GetPointer() ); }

Tensor ExternalClosure::U( const Tensor& /*Lambda*/ ) {
    exit( EXIT_FAILURE );
    return Tensor( 1, 1, 1, -1.0 );
}

void ExternalClosure::DU( Matrix& y, const Vector& Lambda ) { _DU( y.GetPointer(), Lambda.GetPointer() ); }

void ExternalClosure::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    _SolveClosure( lambda.GetPointer(), u.GetPointer(), refLevel );
}

void ExternalClosure::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    _SolveClosure( lambda.GetPointer(), u.GetPointer(), refLevel );
}
