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
        double area = norm( _mesh->GetNormalAtEdge( j ) );
        if( _mesh->BoundaryAtEdge( j ) != NOSLIP && _mesh->BoundaryAtEdge( j ) != DIRICHLET ) {
            fluxFunc( _flux[j], _problem->G( uQ[I], uQ[J], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            fluxFunc( tmp, _problem->G( uQ[J], uQ[I], -_mesh->GetNormalAtEdge( j ) / area, -_mesh->GetNormalAtEdge( j ), level ), level );
            // std::cout << _flux[j] + tmp << std::endl;
            /*
            if( !std::isfinite( ( _flux[j] + tmp )( 0, 0 ) ) ) {
                std::cout << "flux = " << _flux[j] << std::endl;
                std::cout << "uQI " << uQ[I] << std::endl;
                std::cout << "uQJ " << uQ[J] << std::endl;
                exit( EXIT_FAILURE );
            }*/
            // std::cout << "flux = " << _flux[j] << std::endl;
            // std::cout << "uQI " << uQ[I] << std::endl;
            // std::cout << "uQJ " << uQ[J] << std::endl;
            // std::cout << "numCells " << numCells << std::endl;
        }
        else if( _mesh->BoundaryAtEdge( j ) == NOSLIP ) {
            // std::cout << "I = " << I << ", J = " << J << std::endl;
            if( I == numCells ) {
                std::cerr << "ERROR" << std::endl;
                exit( EXIT_FAILURE );
                fluxFunc( _flux[j], _problem->BoundaryFlux( uQ[J], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            }
            else if( J == numCells ) {
                if( I == 9206 ) std::cout << "---> n = " << _mesh->GetNormalAtEdge( j ) << std::endl;
                fluxFunc( _flux[j], _problem->BoundaryFlux( uQ[I], _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            }
            else {
                std::cerr << "ERROR 1" << std::endl;
                exit( EXIT_FAILURE );
            }
        }
    }

    // std::cout << "end edges..." << std::endl;

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
        Matrix rhsFlux( _settings->GetNStates(), _settings->GetNTotalforRefLevel( refLevel[j] ), 0.0 );
        Matrix rhs( _settings->GetNStates(), _settings->GetNTotalforRefLevel( refLevel[j] ), 0.0 );
        auto edges = _mesh->GetEdgesOfCell( j );

        for( unsigned l = 0; l < edges.size(); ++l ) {
            unsigned I = _mesh->CellsAtEdge( edges[l] ).first;
            unsigned J = _mesh->CellsAtEdge( edges[l] ).second;
            // std::cout << "Edge " << edges[l] << " with flux " << _flux[edges[l]] << std::endl;
            // std::cout << _mesh->CellsAtEdge( edges[l] ).first << " " << j << std::endl;
            if( I == cell->GetID() ) {
                // std::cout << "+" << std::endl;
                rhsFlux = rhsFlux + _flux[edges[l]];
                if( neighbors[l] == J && j == 9206 && false ) {
                    Matrix tmp = _flux[j];
                    tmp.reset();
                    fluxFunc( tmp, _problem->G( uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
                    // std::cout << "j == I " << j << " " << I << std::endl;
                    // std::cout << "neighbor is " << neighbors[l] << " " << J << std::endl;
                    std::cout << "diff 1: " << tmp << std::endl;
                    // std::cout << "diff 2: " << _flux[edges[l]] << std::endl;

                    tmp.reset();
                    double area = norm( _mesh->GetNormalAtEdge( j ) );
                    fluxFunc( tmp,
                              _problem->G( uQ[I], uQ[J], _mesh->GetNormalAtEdge( edges[l] ) / area, _mesh->GetNormalAtEdge( edges[l] ), refLevel[j] ),
                              refLevel[j] );
                    std::cout << "diff 2: " << _flux[edges[l]] << std::endl;
                }
            }
            else if( J == cell->GetID() ) {
                // std::cout << "-" << std::endl;1
                rhsFlux = rhsFlux - _flux[edges[l]];
                if( neighbors[l] == I && j == 9206 ) {
                    if( cell->IsBoundaryCell() )
                        std::cout << "Boundary Cell" << std::endl;
                    else
                        std::cout << "No Boundary Cell" << std::endl;
                    if( _mesh->BoundaryAtEdge( j ) == BoundaryType::NOSLIP )
                        std::cout << "Boundary Edge" << std::endl;
                    else
                        std::cout << "No Boundary Edge" << std::endl;
                    Matrix tmp = _flux[j];
                    tmp.reset();
                    // fluxFunc( tmp, _problem->G( uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j]
                    // );
                    fluxFunc( tmp, _problem->BoundaryFlux( uQ[j], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
                    // std::cout << "j == I " << j << " " << I << std::endl;
                    // std::cout << "neighbor is " << neighbors[l] << " " << J << std::endl;
                    std::cout << "diff 1: " << tmp << std::endl;
                    std::cout << "n1: " << cell->GetNormal( l ) << std::endl;
                    // std::cout << "diff 2: " << _flux[edges[l]] << std::endl;

                    tmp.reset();
                    double area = norm( _mesh->GetNormalAtEdge( j ) );
                    fluxFunc( tmp,
                              _problem->G( uQ[I], uQ[J], _mesh->GetNormalAtEdge( edges[l] ) / area, _mesh->GetNormalAtEdge( edges[l] ), refLevel[j] ),
                              refLevel[j] );
                    std::cout << "diff 2: " << _flux[edges[l]] * ( -1.0 ) << std::endl;
                    std::cout << "n2: " << _mesh->GetNormalAtEdge( edges[l] ) << std::endl;
                }
            }
            else {
                std::cerr << "WRONG" << std::endl;
            }
        }

        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells ) {
                fluxFunc( rhs, _problem->BoundaryFlux( uQ[j], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
            else {
                fluxFunc( rhs, _problem->G( uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), refLevel[j] ), refLevel[j] );
            }
        }

        if( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP && j == 9206 ) {
            std::cout << "Cell " << j << std::endl;
            std::cout << rhs - rhsFlux << std::endl;
        }
        // std::cout << rhs - rhsFlux << std::endl;

        // std::cout << std::endl;
        uNew[j] = u[j] - ( dt / cell->GetArea() ) * rhsFlux;
    }
    // exit( EXIT_FAILURE );
    // std::cout << "After advance." << std::endl;
}
