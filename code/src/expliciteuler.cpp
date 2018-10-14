#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Problem* problem ) : TimeSolver( problem ) {}

void ExplicitEuler::Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             std::vector<Matrix>& uNew,
                             std::vector<Matrix>& u,
                             std::vector<Matrix>& uQ ) {
    for( unsigned j = 0; j < _problem->GetMesh()->GetNumCells(); ++j ) {

        auto neighbors = _problem->GetMesh()->GetNeighborIDs( j );
        if( _problem->GetMesh()->GetGrid()[j]->IsBoundaryCell() ) {
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::NOSLIP ) {
                uQ[_problem->GetMesh()->GetNumCells()] = uQ[j];
                for( unsigned k = 0; k < uQ[_problem->GetMesh()->GetNumCells()].columns(); ++k ) {
                    Vector v( 2, 0.0 );
                    v[0] = uQ[_problem->GetMesh()->GetNumCells()]( 1, k );
                    v[1] = uQ[_problem->GetMesh()->GetNumCells()]( 2, k );
                    int index;
                    for( unsigned l = 0; l < neighbors.size(); ++l ) {
                        if( neighbors[l] == _problem->GetMesh()->GetNumCells() ) {
                            index = l;
                        }
                    }
                    Vector n                                       = _problem->GetMesh()->GetUnitNormals( j, index );
                    uQ[_problem->GetMesh()->GetNumCells()]( 1, k ) = -v[0] * n[0];
                    uQ[_problem->GetMesh()->GetNumCells()]( 2, k ) = -v[1] * n[1];
                }
            }
        }
        Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );

        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            rhs += fluxFunc( uQ[j], uQ[neighbors[l]], _problem->GetMesh()->GetUnitNormals( j, l ), _problem->GetMesh()->GetNormals( j, l ) );
        }
        uNew[j] = u[j] - ( _dt / _problem->GetMesh()->GetArea( j ) ) * rhs;
    }
}
