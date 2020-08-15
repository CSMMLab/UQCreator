#include "expliciteuler.h"
#include "mathtools.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh, Problem* problem ) : TimeSolver( settings, mesh, problem ) {
    _cells = _mesh->GetGrid();
    _ghostCell( _settings->GetNStates(), _settings->GetNQuadPoints() );
}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, unsigned )> const& fluxFunc,
                             MatTens& uNew,
                             MatTens& u,
                             MatTens& uQ,
                             double dt,
                             const VectorU& refLevel ) {
    auto numCells = _mesh->GetNumCells();
    Cell* cell;
    VectorU neighbors;
    VectorU edges;
    Matrix uQI( _settings->GetNStates(), _settings->GetNqPE() );
    Matrix uQJ( _settings->GetNStates(), _settings->GetNqPE() );

    // std::cout << "Euler: uQNew = " << uQ[10] << std::endl;

    for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
        std::cout << "ME " << l << std::endl;
// compute flux at edges
#pragma omp parallel for
        for( unsigned j = 0; j < _mesh->GetNEdges(); ++j ) {
            std::cout << "edge = " << j << std::endl;
            _flux[j].reset();    // is this needed? should be fine
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
            std::cout << "ref level is " << level << std::endl;
            std::cout << "nQ at ref level is " << _settings->GetNqPEAtRef( level ) << std::endl;
            for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
                for( unsigned k = 0; k < _settings->GetNqPEAtRef( level ); ++k ) {
                    uQI( s, k ) = uQ[I]( s, l, k );
                    uQJ( s, k ) = uQ[J]( s, l, k );
                }
            }
            std::cout << "written on UQI and UQJ " << std::endl;

            double area = norm( _mesh->GetNormalAtEdge( j ) );
            if( _mesh->BoundaryAtEdge( j ) != NOSLIP && _mesh->BoundaryAtEdge( j ) != DIRICHLET ) {
                std::cout << "volume " << std::endl;
                fluxFunc( _flux[j], _problem->G( uQI, uQJ, _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
            }
            else if( _mesh->BoundaryAtEdge( j ) == NOSLIP ) {
                std::cout << "noslip " << std::endl;
                if( I == numCells ) {
                    std::cerr << "ERROR" << std::endl;
                    exit( EXIT_FAILURE );
                    fluxFunc(
                        _flux[j], _problem->BoundaryFlux( uQJ, _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
                }
                else if( J == numCells ) {
                    std::cout << j << std::endl;
                    fluxFunc(
                        _flux[j], _problem->BoundaryFlux( uQI, _mesh->GetNormalAtEdge( j ) / area, _mesh->GetNormalAtEdge( j ), level ), level );
                    std::cout << "Flux function computed." << std::endl;
                }
                else {
                    std::cerr << "ERROR 1" << std::endl;
                    exit( EXIT_FAILURE );
                }
            }
        }
        std::cout << "flux computation DONE " << std::endl;

#pragma omp parallel for
        for( unsigned j = 0; j < numCells; ++j ) {
            cell      = _cells[j];
            neighbors = cell->GetNeighborIDs();    // neighbors at cell j

            if( cell->IsBoundaryCell() ) {
                if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                    uNew[j] = uQ[j];
                    continue;
                }
            }
            Matrix rhs( _settings->GetNStates(), _settings->GetNqPEAtRef( refLevel[j] ), 0.0 );
            edges = _mesh->GetEdgesOfCell( j );

            for( unsigned ngh = 0; ngh < edges.size(); ++ngh ) {
                unsigned I = _mesh->CellsAtEdge( edges[ngh] ).first;
                unsigned J = _mesh->CellsAtEdge( edges[ngh] ).second;
                if( I == cell->GetID() ) {
                    rhs = rhs + _flux[edges[ngh]];
                }
                else if( J == cell->GetID() ) {
                    rhs = rhs - _flux[edges[ngh]];
                }
                else {
                    std::cerr << "WRONG" << std::endl;
                    std::cout << "Edges of Cell " << j << std::endl;
                    std::cout << "Cell index is " << cell->GetID() << std::endl;
                    std::cout << I << " " << J << std::endl;
                }
            }

            // TODO: Check if refLevel and uQ[j] is correct
            for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
                for( unsigned k = 0; k < _settings->GetNqPEAtRef( refLevel[j] ); ++k ) {
                    uQJ( s, k ) = uQ[j]( s, l, k );
                }
            }

            auto uQNew = uQJ - ( dt / cell->GetArea() ) * rhs;

            for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
                for( unsigned k = 0; k < uQNew.columns(); ++k ) {
                    uNew[j]( s, l, k ) = uQNew( s, k );
                }
            }
        }
        std::cout << "up uQ DONE " << std::endl;
    }
}
