#include "problems/externalproblem.h"

ExternalProblem::ExternalProblem( double* ( *G )(double*, double*, double*, double*, unsigned),
                                  double* ( *F )(double*),
                                  double ( *ComputeDt )( unsigned ),
                                  double* ( *IC )(double*, double*),
                                  Settings* settings )
    : Problem( settings ), _G( G ), _F( F ), _ComputeDt( ComputeDt ), _IC( IC ) {}

ExternalProblem::~ExternalProblem() {}

void ExternalProblem::Solve() {}

Vector ExternalProblem::G( const Vector& /*u*/, const Vector& /*v*/, const Vector& /*nUnit*/, const Vector& /*n*/ ) {
    exit( EXIT_FAILURE );
    return Vector( 1, -1.0 );
}

Matrix ExternalProblem::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    double* res = _G( u.GetPointer(), v.GetPointer(), nUnit.GetPointer(), n.GetPointer(), level );
    return Matrix( _settings->GetNStates(), _settings->GetNQuadPoints(), res );
}

Matrix ExternalProblem::F( const Vector& u ) const {
    double* res = _F( u.GetPointer() );
    return Matrix( _settings->GetNStates(), _settings->GetMeshDimension(), res );
}

Matrix ExternalProblem::F( const Matrix& /*u*/ ) {
    exit( EXIT_FAILURE );
    return Matrix( 1, 1, -1.0 );
}

Matrix ExternalProblem::BoundaryFlux( const Matrix& /*u*/, const Vector& /*nUnit*/, const Vector& /*n*/, unsigned /*level*/ ) const {
    exit( EXIT_FAILURE );
    return Matrix( 1, 1, -1.0 );
}

double ExternalProblem::ComputeDt( const Tensor& /*u*/, double dx, unsigned /*level*/ ) const { return _ComputeDt( dx ); }

Vector ExternalProblem::IC( const Vector& x, const Vector& xi ) {
    double* res = _IC( x.GetPointer(), xi.GetPointer() );
    return Vector( _settings->GetNStates(), res );
}

Vector ExternalProblem::LoadIC( const Vector& /*x*/, const Vector& /*xi*/ ) {
    exit( EXIT_FAILURE );
    return Vector( 1, -1.0 );
}
