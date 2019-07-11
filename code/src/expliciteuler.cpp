#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh, Problem* problem ) : TimeSolver( settings, mesh, problem ) {
    _cells = _mesh->GetGrid();
    _ghostCell( _settings->GetNStates(), _settings->GetNQuadPoints() );
}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, unsigned )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& u,
                             MatVec& uQ,
                             double dt,
                             const VectorU& refLevel ) {
    auto numCells = _mesh->GetNumCells();

#pragma omp parallel for
    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell     = _cells[j];
        auto neighbors = cell->GetNeighborIDs();
        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
        }

        Matrix rhs( _settings->GetNStates(), _settings->GetNTotalforRefLevel( refLevel[j] ), 0.0 );
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells ) {
                fluxFunc( rhs, _problem->BoundaryFlux( uQ[j], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
            else {
                fluxFunc( rhs, _problem->G( uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
        }
        uNew[j] = u[j] - ( dt / cell->GetArea() ) * rhs;
    }
}
