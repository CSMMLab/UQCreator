#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& uQ,
                             double dt ) {
    auto numCells = _mesh->GetNumCells();
    auto cells    = _mesh->GetGrid();
    Matrix ghostCell( _settings->GetNStates(), _settings->GetNQuadPoints() );

#pragma omp parallel for private( ghostCell )
    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell     = cells[j];
        auto neighbors = cell->GetNeighborIDs();
        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = uQ[j];
                continue;
            }
            else if( cell->GetBoundaryType() == BoundaryType::NOSLIP ) {
                ghostCell = uQ[j];
                Vector v( 2, 0.0 );
                v.reset();
                v[0]              = ghostCell( 1, 0 ) / ghostCell( 0, 0 );
                v[1]              = ghostCell( 2, 0 ) / ghostCell( 0, 0 );
                Vector n          = cell->GetBoundaryUnitNormal();
                double vn         = dot( n, v );
                Vector Vn         = vn * n;
                Vector Vb         = -Vn + v;
                double velMagB    = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                double velMag     = v[0] * v[0] + v[1] * v[1];
                double rho        = ghostCell( 0, 0 );
                ghostCell( 1, 0 ) = rho * ( Vb[0] );
                ghostCell( 2, 0 ) = rho * ( Vb[1] );
                ghostCell( 3, 0 ) += rho * 0.5 * ( velMagB - velMag );
            }
            else if( cell->GetBoundaryType() == BoundaryType::SWWALL ) {
                ghostCell = uQ[j];
                Vector v( 2, 0.0 );
                v[0]              = ghostCell( 1, 0 ) / ghostCell( 0, 0 );
                v[1]              = ghostCell( 2, 0 ) / ghostCell( 0, 0 );
                Vector n          = cell->GetBoundaryUnitNormal();
                double vn         = dot( n, v );
                Vector Vn         = vn * n;
                Vector Vb         = -Vn + v;
                double rho        = ghostCell( 0, 0 );
                ghostCell( 1, 0 ) = rho * ( Vb[0] );
                ghostCell( 2, 0 ) = rho * ( Vb[1] );
            }
        }
        bool boundary = false;
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells )
                boundary = true;
        }
        Matrix rhs( uQ[0].rows(), uQ[0].columns(), 0.0 );
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells ) {
                fluxFunc( rhs, ghostCell, ghostCell, cell->GetUnitNormal( l ), cell->GetNormal( l ) );
            }
            else {
                fluxFunc( rhs, uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ) );
            }
        }
        uNew[j] = uQ[j] - ( dt / cell->GetArea() ) * rhs;
    }
}
