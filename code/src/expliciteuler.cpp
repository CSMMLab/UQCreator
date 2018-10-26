#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Settings* settings, Mesh* mesh ) : TimeSolver( settings, mesh ) {}

void ExplicitEuler::Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             std::vector<Matrix>& uNew,
                             std::vector<Matrix>& u,
                             std::vector<Matrix>& uQ ) {
    for( unsigned j = 0; j < _mesh->GetNumCells(); ++j ) {

        auto neighbors = _mesh->GetNeighborIDs( j );
        if( _mesh->GetGrid()[j]->IsBoundaryCell() ) {
            if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            if( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP ) {
                uQ[_mesh->GetNumCells()] = uQ[j];
                for( unsigned k = 0; k < uQ[_mesh->GetNumCells()].columns(); ++k ) {
                    Vector v( 2, 0.0 );
                    v[0]           = uQ[_mesh->GetNumCells()]( 1, k ) / uQ[_mesh->GetNumCells()]( 0, k );
                    v[1]           = uQ[_mesh->GetNumCells()]( 2, k ) / uQ[_mesh->GetNumCells()]( 0, k );
                    unsigned index = 100;
                    for( unsigned l = 0; l < neighbors.size(); ++l ) {
                        if( neighbors[l] == _mesh->GetNumCells() ) {
                            index = l;
                        }
                    }
                    if( index == 100 ) {
                        std::cerr << "Boundary Cell " << j << " has no ghost cell neighbor." << std::endl;
                        exit( EXIT_FAILURE );
                    }
                    Vector n                         = _mesh->GetUnitNormals( j, index );
                    double vn                        = n[0] * v[0] + n[1] * v[1];
                    Vector Vn                        = vn * n;
                    Vector Vb                        = -Vn + v;
                    double velMagB                   = Vb[0] * Vb[0] + Vb[1] * Vb[1];
                    double velMag                    = v[0] * v[0] + v[1] * v[1];
                    double rho                       = uQ[_mesh->GetNumCells()]( 0, k );
                    uQ[_mesh->GetNumCells()]( 1, k ) = rho * ( Vb[0] );
                    uQ[_mesh->GetNumCells()]( 2, k ) = rho * ( Vb[1] );
                    uQ[_mesh->GetNumCells()]( 3, k ) += rho * 0.5 * ( velMagB - velMag );
                }
            }
        }
        Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );

        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            if( _mesh->GetBoundaryType( j ) == BoundaryType::NOSLIP && neighbors[l] == _mesh->GetNumCells() ) {
                rhs += fluxFunc( uQ[neighbors[l]], uQ[neighbors[l]], _mesh->GetUnitNormals( j, l ), _mesh->GetNormals( j, l ) );
            }
            else {
                rhs += fluxFunc( uQ[j], uQ[neighbors[l]], _mesh->GetUnitNormals( j, l ), _mesh->GetNormals( j, l ) );
            }
        }
        uNew[j] = u[j] - ( _dt / _mesh->GetArea( j ) ) * rhs;
    }
}
