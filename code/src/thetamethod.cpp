#include "thetamethod.h"

ThetaMethod::ThetaMethod( Problem* problem, double theta ) : TimeSolver( problem ), _theta( theta ) {}

void ThetaMethod::Advance(
    std::function<Matrix( const Matrix&, const Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
    std::vector<Matrix>& uNew,
    std::vector<Matrix>& u,
    std::vector<Matrix>& lambda ) {
    if( _theta == 0 ) {    // explicit Euler
        for( unsigned j = 0; j < _problem->GetMesh()->GetNumCells(); ++j ) {

            blaze::DynamicVector<unsigned> neighbors = _problem->GetMesh()->GetNeighborsIndex( j );
            if( _problem->GetMesh()->GetGrid()[j]->IsBoundaryCell() ) {
                lambda[_problem->GetMesh()->GetNumCells()] = lambda[j];
                std::cout << "Boundary Cell " << j << ", with values " << lambda[j] << std::endl;
                std::cout << neighbors[0] << " " << lambda[neighbors[0]] << std::endl;
                std::cout << neighbors[1] << " " << lambda[neighbors[1]] << std::endl;
            }
            Matrix rhs( lambda[0].rows(), lambda[0].columns(), 0.0 );

            for( unsigned l = 0; l < neighbors.size(); ++l ) {
                rhs += fluxFunc( lambda[j],
                                 lambda[j],
                                 lambda[neighbors[l]],
                                 lambda[j],
                                 _problem->GetMesh()->GetUnitNormals( j, l ),
                                 _problem->GetMesh()->GetNormals( j, l ) );
            }
            uNew[j] = u[j] - ( _dt / _problem->GetMesh()->GetArea( j ) ) * rhs;
        }
    }
    else if( _theta == 0.5 ) {    // Crank-Nicolson
    }
    else if( _theta == 1 ) {    // implicit Euler
    }
}
