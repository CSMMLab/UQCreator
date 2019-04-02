#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {
    _cells = _mesh->GetGrid();
    _ghostCell( _settings->GetNStates(), _settings->GetNQuadPoints() );
}

void ExplicitEuler::Advance( std::function<void( Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector&, unsigned )> const& fluxFunc,
                             MatVec& uNew,
                             MatVec& u,
                             MatVec& uQ,
                             double dt,
                             const VectorU& refLevel ) {
    auto numCells = _mesh->GetNumCells();

#pragma omp parallel for private( _ghostCell )
    for( unsigned j = 0; j < numCells; ++j ) {
        Cell* cell     = _cells[j];
        auto neighbors = cell->GetNeighborIDs();
        if( cell->IsBoundaryCell() ) {
            if( cell->GetBoundaryType() == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            else if( cell->GetBoundaryType() == BoundaryType::NOSLIP ) {
                _ghostCell = uQ[j];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < uQ[numCells].columns(); ++k ) {
                    v.reset();
                    v[0]               = _ghostCell( 1, k ) / _ghostCell( 0, k );
                    v[1]               = _ghostCell( 2, k ) / _ghostCell( 0, k );
                    Vector n           = cell->GetBoundaryUnitNormal();
                    double vn          = dot( n, v );
                    Vector Vn          = vn * n;
                    Vector Vb          = -Vn + v;
                    double velMagB     = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                    double velMag      = v[0] * v[0] + v[1] * v[1];
                    double rho         = _ghostCell( 0, k );
                    _ghostCell( 1, k ) = rho * ( Vb[0] );
                    _ghostCell( 2, k ) = rho * ( Vb[1] );
                    _ghostCell( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                }
            }
            else if( cell->GetBoundaryType() == BoundaryType::SWWALL ) {
                _ghostCell = uQ[j];
                Vector v( 2, 0.0 );
                for( unsigned k = 0; k < uQ[numCells].columns(); ++k ) {
                    v[0]               = _ghostCell( 1, k ) / _ghostCell( 0, k );
                    v[1]               = _ghostCell( 2, k ) / _ghostCell( 0, k );
                    Vector n           = cell->GetBoundaryUnitNormal();
                    double vn          = dot( n, v );
                    Vector Vn          = vn * n;
                    Vector Vb          = -Vn + v;
                    double rho         = _ghostCell( 0, k );
                    _ghostCell( 1, k ) = rho * ( Vb[0] );
                    _ghostCell( 2, k ) = rho * ( Vb[1] );
                }
            }
        }

        Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );
        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( ( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP || _mesh->GetBoundaryType( j ) == BoundaryType::SWWALL ) &&
                neighbors[l] == numCells ) {
                fluxFunc(
                    rhs, _ghostCell, _ghostCell, cell->GetUnitNormal( l ), cell->GetNormal( l ), _settings->GetNTotalforRefLevel( refLevel[j] ) );
            }
            else {
                fluxFunc(
                    rhs, uQ[j], uQ[neighbors[l]], cell->GetUnitNormal( l ), cell->GetNormal( l ), _settings->GetNTotalforRefLevel( refLevel[j] ) );
            }
        }
        uNew[j] = u[j] - ( dt / cell->GetArea() ) * rhs;
    }
}
