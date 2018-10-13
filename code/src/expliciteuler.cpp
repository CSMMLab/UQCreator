#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Problem* problem ) : TimeSolver( problem ) {}

void ExplicitEuler::Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             std::vector<Matrix>& uNew,
                             std::vector<Matrix>& u,
                             std::vector<Matrix>& lambda ) {
    for( unsigned j = 0; j < _problem->GetMesh()->GetNumCells(); ++j ) {

        blaze::DynamicVector<unsigned> neighbors = _problem->GetMesh()->GetNeighborsIndex( j );
        if( _problem->GetMesh()->GetGrid()[j]->IsBoundaryCell() ) {
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::NEUMANN ) {
                lambda[_problem->GetMesh()->GetNumCells()] = lambda[j];
                // lambda[_problem->GetMesh()->GetNumCells()](1,:) *= -1.0;
            }
        }
        Matrix rhs( lambda[0].rows(), lambda[0].columns(), 0.0 );

        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            rhs += fluxFunc( lambda[j], lambda[neighbors[l]], _problem->GetMesh()->GetUnitNormals( j, l ), _problem->GetMesh()->GetNormals( j, l ) );
        }
        uNew[j] = u[j] - ( _dt / _problem->GetMesh()->GetArea( j ) ) * rhs;
    }
}
