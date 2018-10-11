#include "thetamethod.h"

ThetaMethod::ThetaMethod( Problem* problem, double theta ) : TimeSolver( problem ), _theta( theta ) {}

void ThetaMethod::Advance(
    std::function<Matrix( const Matrix&, const Matrix&, const Matrix&, const Matrix&, const Vector&, const Vector& )> const& fluxFunc,
    std::vector<Matrix>& uNew,
    std::vector<Matrix>& u,
    std::vector<Matrix>& lambda ) {
    if( _theta == 0 ) {    // explicit Euler
        for( unsigned j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
            blaze::DynamicVector<unsigned> neighbors = _problem->GetMesh()->GetNeighbors( j );
            Matrix rhs( lambda[0].rows(), lambda[0].columns(), 0.0 );

            for( unsigned l = 0; l < neighbors.size(); ++l ) {
                rhs += fluxFunc( lambda[j - 1],
                                 lambda[j],
                                 lambda[neighbors[l]],
                                 lambda[j + 2],
                                 _problem->GetMesh()->GetUnitNormals( j, l ),
                                 _problem->GetMesh()->GetNormals( j, l ) );
            }
            uNew[j] = u[j] - _dt * rhs;
        }
    }
    else if( _theta == 0.5 ) {    // Crank-Nicolson
    }
    else if( _theta == 1 ) {    // implicit Euler
    }
}
