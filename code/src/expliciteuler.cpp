#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& u,
                             MatVec& uQ ) {
    auto numCells = _mesh->GetNumCells();
    auto edges    = _mesh->GetEdges();
    auto cells    = _mesh->GetGrid();
    std::vector<Matrix> rhs( numCells, Matrix( _settings->GetNStates(), _settings->GetNQuadPoints() ) );

    //#pragma omp parallel for private( ghostCell ) reduction(+:rhs)
    for( unsigned i = 0; i < edges.size(); ++i ) {
        Cell* cellA      = edges[i].first;
        Cell* cellB      = edges[i].second;
        auto neighborIDs = cellA->GetNeighborIDs();

        if( cellB != nullptr ) {
            unsigned faceID = 0;
            for( unsigned j = 0; j < neighborIDs.size(); ++j ) {
                if( neighborIDs[j] == cellB->GetID() ) {
                    faceID = j;
                    break;
                }
            }
            assert( neighborIDs[faceID] == cellB->GetID() );
            if( cellA->GetBoundaryType() == BoundaryType::NOSLIP ) {
                Matrix ghostCell = uQ[cellA->GetID()];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < _settings->GetNMoments(); ++k ) {
                    v.reset();
                    v[0]              = ghostCell( 1, k ) / ghostCell( 0, k );
                    v[1]              = ghostCell( 2, k ) / ghostCell( 0, k );
                    Vector n          = cellA->GetBoundaryUnitNormal();
                    double vn         = dot( n, v );
                    Vector Vn         = vn * n;
                    Vector Vb         = -Vn + v;
                    double velMagB    = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                    double velMag     = v[0] * v[0] + v[1] * v[1];
                    double rho        = ghostCell( 0, k );
                    ghostCell( 1, k ) = rho * ( Vb[0] );
                    ghostCell( 2, k ) = rho * ( Vb[1] );
                    ghostCell( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                }
                Matrix flux( _settings->GetNStates(), _settings->GetNQuadPoints() );
                fluxFunc( flux, ghostCell, ghostCell, cellA->GetUnitNormal( faceID ), cellA->GetNormal( faceID ) );
                rhs[cellA->GetID()] += flux;
                continue;
            }
            else if( cellB->GetBoundaryType() == BoundaryType::NOSLIP ) {
                Matrix ghostCell = uQ[cellB->GetID()];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < _settings->GetNMoments(); ++k ) {
                    v.reset();
                    v[0]              = ghostCell( 1, k ) / ghostCell( 0, k );
                    v[1]              = ghostCell( 2, k ) / ghostCell( 0, k );
                    Vector n          = cellB->GetBoundaryUnitNormal();
                    double vn         = dot( n, v );
                    Vector Vn         = vn * n;
                    Vector Vb         = -Vn + v;
                    double velMagB    = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                    double velMag     = v[0] * v[0] + v[1] * v[1];
                    double rho        = ghostCell( 0, k );
                    ghostCell( 1, k ) = rho * ( Vb[0] );
                    ghostCell( 2, k ) = rho * ( Vb[1] );
                    ghostCell( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                }
                Matrix flux( _settings->GetNStates(), _settings->GetNQuadPoints() );
                fluxFunc( flux, ghostCell, ghostCell, cellB->GetUnitNormal( faceID ), cellB->GetNormal( faceID ) );
                rhs[cellB->GetID()] += flux;
                continue;
            }

            Matrix flux( _settings->GetNStates(), _settings->GetNQuadPoints() );
            fluxFunc( flux, uQ[cellA->GetID()], uQ[cellB->GetID()], cellA->GetUnitNormal( faceID ), cellA->GetNormal( faceID ) );
            rhs[cellA->GetID()] += flux;
            rhs[cellB->GetID()] -= flux;
        }
        else {
            unsigned faceID = 0;
            for( unsigned j = 0; j < neighborIDs.size(); ++j ) {
                if( neighborIDs[j] == numCells ) {
                    faceID = j;
                    break;
                }
            }
            assert( neighborIDs[faceID] == numCells );
            if( cellA->GetBoundaryType() == BoundaryType::NOSLIP ) {
                Matrix ghostCell = uQ[cellA->GetID()];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < _settings->GetNMoments(); ++k ) {
                    v.reset();
                    v[0]              = ghostCell( 1, k ) / ghostCell( 0, k );
                    v[1]              = ghostCell( 2, k ) / ghostCell( 0, k );
                    Vector n          = cellA->GetBoundaryUnitNormal();
                    double vn         = dot( n, v );
                    Vector Vn         = vn * n;
                    Vector Vb         = -Vn + v;
                    double velMagB    = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                    double velMag     = v[0] * v[0] + v[1] * v[1];
                    double rho        = ghostCell( 0, k );
                    ghostCell( 1, k ) = rho * ( Vb[0] );
                    ghostCell( 2, k ) = rho * ( Vb[1] );
                    ghostCell( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                }
                Matrix flux( _settings->GetNStates(), _settings->GetNQuadPoints() );
                fluxFunc( flux, ghostCell, ghostCell, cellA->GetUnitNormal( faceID ), cellA->GetNormal( faceID ) );
                rhs[cellA->GetID()] += flux;
            }
        }
    }

    //#pragma omp parallel for
    for( unsigned j = 0; j < numCells; ++j ) {
        uNew[j] = u[j] - ( _dt / cells[j]->GetArea() ) * rhs[j];
    }

    /*
    #pragma omp parallel for private( ghostCell )
        for( unsigned j = 0; j < numCells; ++j ) {
            Cell* cell     = cells[j];
            auto neighbors = cell->GetNeighborIDs();
            if( cell->IsBoundaryCell() ) {
                if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                    uNew[j] = u[j];
                    continue;
                }
                else if( cell->GetBoundaryType() == BoundaryType::NOSLIP ) {
                    ghostCell = uQ[j];
                    Vector v( 2, 0.0 );
                    for( unsigned k = 0; k < uQ[numCells].columns(); ++k ) {
                        v.reset();
                        v[0]              = ghostCell( 1, k ) / ghostCell( 0, k );
                        v[1]              = ghostCell( 2, k ) / ghostCell( 0, k );
                        Vector n          = cell->GetBoundaryUnitNormal();
                        double vn         = dot( n, v );
                        Vector Vn         = vn * n;
                        Vector Vb         = -Vn + v;
                        double velMagB    = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                        double velMag     = v[0] * v[0] + v[1] * v[1];
                        double rho        = ghostCell( 0, k );
                        ghostCell( 1, k ) = rho * ( Vb[0] );
                        ghostCell( 2, k ) = rho * ( Vb[1] );
                        ghostCell( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                    }
                }
            }

            Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );
            for( unsigned l = 0; l < neighbors.size(); ++l ) {
                // std::cout << cell->GetID() << "\t" << neighbors[l] << std::endl;
                if( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP && neighbors[l] == numCells ) {
                    fluxFunc( rhs, ghostCell, ghostCell, cell->GetUnitNormal( l ), cell->GetNormal( l ) );
                }
                else {
                    fluxFunc( rhs, uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ) );
                }
            }
            uNew[j] = u[j] - ( _dt / cell->GetArea() ) * rhs;
        }
        */
}
