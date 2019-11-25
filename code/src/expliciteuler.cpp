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
    // std::cout << "In advance..." << std::endl;
// compute flux at edges
#pragma omp parallel for
    for( unsigned j = 0; j < _mesh->GetNEdges(); ++j ) {
        _flux[j].reset();    // is this needed? should be fine
        Matrix tmp      = _flux[j];
        std::pair cells = _mesh->CellsAtEdge( j );
        unsigned I      = cells.first;
        unsigned J      = cells.second;
        // std::cout << "edge " << j << ": " << I << " " << J << std::endl;
        unsigned level = refLevel[I];    // take max ref level of cells I and J
        if( level < refLevel[J] ) level = refLevel[J];
        double area = norm( _mesh->GetNormalsAtEdge( j ) );
        if( _mesh->BoundaryAtEdge( j ) != NOSLIP && _mesh->BoundaryAtEdge( j ) != DIRICHLET ) {
            fluxFunc( _flux[j], _problem->G( uQ[I], uQ[J], _mesh->GetNormalsAtEdge( j ) / area, _mesh->GetNormalsAtEdge( j ), level ), level );
            fluxFunc( tmp, _problem->G( uQ[J], uQ[I], -_mesh->GetNormalsAtEdge( j ) / area, -_mesh->GetNormalsAtEdge( j ), level ), level );
            std::cout << _flux[j] + tmp << std::endl;
            if( !std::isfinite( ( _flux[j] + tmp )( 0, 0 ) ) ) {
                std::cout << "flux = " << _flux[j] << std::endl;
                std::cout << "uQI " << uQ[I] << std::endl;
                std::cout << "uQJ " << uQ[J] << std::endl;
                exit( EXIT_FAILURE );
            }
            // std::cout << "flux = " << _flux[j] << std::endl;
            // std::cout << "uQI " << uQ[I] << std::endl;
            // std::cout << "uQJ " << uQ[J] << std::endl;
            // std::cout << "numCells " << numCells << std::endl;
        }
        else if( _mesh->BoundaryAtEdge( j ) == NOSLIP ) {
            fluxFunc( _flux[j], _problem->BoundaryFlux( uQ[I], _mesh->GetNormalsAtEdge( j ) / area, _mesh->GetNormalsAtEdge( j ), level ), level );
        }
    }

    // std::cout << "end edges..." << std::endl;

#pragma omp parallel for
    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell = _cells[j];

        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
        }
        Matrix rhs( _settings->GetNStates(), _settings->GetNTotalforRefLevel( refLevel[j] ), 0.0 );
        /*

                auto edges = _mesh->GetEdgesOfCell( j );
                // std::cout << "Cell " << j << ": ";
                for( unsigned l = 0; l < edges.size(); ++l ) {
                    // std::cout << "Edge " << edges[l] << " with flux " << _flux[edges[l]] << std::endl;
                    // std::cout << _mesh->CellsAtEdge( edges[l] ).first << " " << j << std::endl;
                    if( _mesh->CellsAtEdge( edges[l] ).first == cell->GetID() ) {
                        // std::cout << "+" << std::endl;
                        rhs = rhs + _flux[edges[l]];
                    }
                    else {
                        // std::cout << "-" << std::endl;1
                        rhs = rhs - _flux[edges[l]];
                    }
                }*/

        auto neighbors = cell->GetNeighborIDs();

        rhs.reset();
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells ) {
                fluxFunc( rhs, _problem->BoundaryFlux( uQ[j], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
            else {
                fluxFunc( rhs, _problem->G( uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
        }

        // std::cout << std::endl;
        uNew[j] = u[j] - ( dt / cell->GetArea() ) * rhs;
    }
    // exit( EXIT_FAILURE );
    // std::cout << "After advance." << std::endl;
}
