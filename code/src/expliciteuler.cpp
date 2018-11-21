#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {}

void ExplicitEuler::Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& u,
                             MatVec& uQ ) {
    auto numCells = _mesh->GetNumCells();
    auto edges    = _mesh->GetEdges();
    auto cells    = _mesh->GetGrid();
    Matrix ghostCell( _settings->GetNStates(), _settings->GetNQuadPoints() );
    FluxMatrix fluxes( numCells );
    std::vector<Matrix> rhs( numCells, Matrix( _settings->GetNStates(), _settings->GetNMoments() ) );
    for( const auto& edge : edges ) {
        unsigned l     = 0;
        auto neighbors = edge.first->GetNeighborIDs();
        for( ; l < neighbors.size(); ++l ) {
            if( neighbors[l] == edge.second->GetID() ) break;
        }
        auto flux = fluxFunc( uQ[edge.first->GetID()], uQ[edge.second->GetID()], edge.first->GetUnitNormal( l ), edge.first->GetNormal( l ) );
        rhs[edge.first->GetID()] += flux;
        rhs[edge.second->GetID()] += -flux;
    }

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
