#include "expliciteuler.h"

ExplicitEuler::ExplicitEuler( Problem* problem ) : TimeSolver( problem ) {}

void ExplicitEuler::Advance( std::function<Matrix( const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
                             std::vector<Matrix>& uNew,
                             std::vector<Matrix>& u,
                             std::vector<Matrix>& uQ ) {
    for( unsigned j = 0; j < _problem->GetMesh()->GetNumCells(); ++j ) {

        auto neighbors = _problem->GetMesh()->GetNeighborIDs( j );
        if( _problem->GetMesh()->GetGrid()[j]->IsBoundaryCell() ) {
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                uNew[j] = u[j];
                continue;
            }
            if( _problem->GetMesh()->GetBoundaryType( j ) == BoundaryType::NOSLIP ) {
                uQ[_problem->GetMesh()->GetNumCells()] = uQ[j];
                for( unsigned k = 0; k < uQ[_problem->GetMesh()->GetNumCells()].columns(); ++k ) {
                    Vector v( 2, 0.0 );
                    v[0]           = uQ[_problem->GetMesh()->GetNumCells()]( 1, k ) / uQ[_problem->GetMesh()->GetNumCells()]( 0, k );
                    v[1]           = uQ[_problem->GetMesh()->GetNumCells()]( 2, k ) / uQ[_problem->GetMesh()->GetNumCells()]( 0, k );
                    unsigned index = 100;
                    // if( _problem->GetMesh()->GetGrid()[j]->IsBoundaryCell() ) {
                    //    std::cout << "Is boundary cell and has " << neighbors.size() << " neighbors" << std::endl;
                    //}
                    for( unsigned l = 0; l < neighbors.size(); ++l ) {
                        // std::cout << neighbors[l] << " != " << _problem->GetMesh()->GetNumCells() << std::endl;
                        if( neighbors[l] == _problem->GetMesh()->GetNumCells() ) {
                            index = l;
                        }
                    }
                    if( index == 100 ) {
                        std::cerr << "Boundary Cell " << j << " has no ghost cell neighbor." << std::endl;
                        exit( EXIT_FAILURE );
                    }
                    Vector n                                       = _problem->GetMesh()->GetUnitNormals( j, index );
                    double vn                                      = n[0] * v[0] + n[1] * v[1];
                    Vector Vn                                      = vn * n;
                    Vector Vt                                      = v - Vn;
                    uQ[_problem->GetMesh()->GetNumCells()]( 1, k ) = uQ[_problem->GetMesh()->GetNumCells()]( 0, k ) * ( -Vn[0] + Vt[0] );
                    uQ[_problem->GetMesh()->GetNumCells()]( 2, k ) = uQ[_problem->GetMesh()->GetNumCells()]( 0, k ) * ( -Vn[1] + Vt[1] );
                }
            }
        }
        Matrix rhs( u[0].rows(), u[0].columns(), 0.0 );

        for( unsigned l = 0; l < neighbors.size(); ++l ) {
            rhs += fluxFunc( uQ[j], uQ[neighbors[l]], _problem->GetMesh()->GetUnitNormals( j, l ), _problem->GetMesh()->GetNormals( j, l ) );
        }
        uNew[j] = u[j] - ( _dt / _problem->GetMesh()->GetArea( j ) ) * rhs;
    }
}
