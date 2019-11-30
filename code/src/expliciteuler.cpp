#include "expliciteuler.h"
#include "mathtools.h"

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
// compute flux at edges
#pragma omp parallel for
    for( unsigned j = 0; j < _mesh->GetNEdges(); ++j ) {
        _flux[j].reset();    // is this needed? should be fine
        Matrix tmp      = _flux[j];
        std::pair cells = _mesh->CellsAtEdge( j );
        unsigned I      = cells.first;
        unsigned J      = cells.second;
        unsigned level;

        if( I == numCells )    // ref level at numCells does not exist
            level = refLevel[J];
        else if( J == numCells )
            level = refLevel[I];
        else {
            level = refLevel[I];    // take max ref level of cells I and J
            if( level < refLevel[J] ) level = refLevel[J];
        }
        double area = norm( _mesh->GetNormalAtEdge( j ) );
        if( _mesh->BoundaryAtEdge( j ) != NOSLIP && _mesh->BoundaryAtEdge( j ) != DIRICHLET ) {
            fluxFunc( _flux[j], _problem->G( uQ[I], uQ[J], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            fluxFunc( tmp, _problem->G( uQ[J], uQ[I], -_mesh->GetNormalAtEdge( j ) / area, -_mesh->GetNormalAtEdge( j ), level ), level );
        }
        else if( _mesh->BoundaryAtEdge( j ) == NOSLIP ) {
            // std::cout << "I = " << I << ", J = " << J << std::endl;
            if( I == numCells ) {
                std::cerr << "ERROR" << std::endl;
                exit( EXIT_FAILURE );
                fluxFunc( _flux[j], _problem->BoundaryFlux( uQ[J], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            }
            else if( J == numCells ) {
                fluxFunc( _flux[j], _problem->BoundaryFlux( uQ[I], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            }
            else {
                std::cerr << "ERROR 1" << std::endl;
                exit( EXIT_FAILURE );
            }
        }
    }

#pragma omp parallel for
    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell     = _cells[j];
        auto neighbors = cell->GetNeighborIDs();    // neighbors at cell j

        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
        }
        Matrix rhs( _settings->GetNStates(), _settings->GetNTotalforRefLevel( refLevel[j] ), 0.0 );
        auto edges = _mesh->GetEdgesOfCell( j );

        for( unsigned l = 0; l < edges.size(); ++l ) {
            unsigned I = _mesh->CellsAtEdge( edges[l] ).first;
            unsigned J = _mesh->CellsAtEdge( edges[l] ).second;
            if( I == cell->GetID() ) {
                rhs = rhs + _flux[edges[l]];
            }
            else if( J == cell->GetID() ) {
                rhs = rhs - _flux[edges[l]];
            }
            else {
                std::cerr << "WRONG" << std::endl;
            }
        }

        uNew[j] = u[j] - ( dt / cell->GetArea() ) * rhs;
    }
}
