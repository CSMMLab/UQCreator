#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& u,
                             MatVec& uQ ) {
    auto numCells = _mesh->GetNumCells();
    auto cells    = _mesh->GetGrid();
    Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );

    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell     = cells[j];
        auto neighbors = cell->GetNeighborIDs();
        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            else if( cell->GetBoundaryType() == BoundaryType::NOSLIP ) {
                uQ[numCells] = uQ[j];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < uQ[numCells].columns(); ++k ) {
                    v.reset();
                    v[0]      = uQ[numCells]( 1, k ) / uQ[numCells]( 0, k );
                    v[1]      = uQ[numCells]( 2, k ) / uQ[numCells]( 0, k );
                    Vector n  = cell->GetBoundaryUnitNormal();
                    double vn = dot( n, v );
                    Vector Vn = vn * n;
                    Vector Vt = v - Vn;

                    uQ[numCells]( 1, k ) = uQ[numCells]( 0, k ) * ( -Vn[0] + Vt[0] );
                    uQ[numCells]( 2, k ) = uQ[numCells]( 0, k ) * ( -Vn[1] + Vt[1] );
                }
            }
        }

        rhs.reset();
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            fluxFunc( rhs, uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ) );
        }
        uNew[j] = u[j] - ( _dt / cell->GetArea() ) * rhs;
    }
}
