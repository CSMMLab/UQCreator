#include "thetamethod.h"

ThetaMethod::ThetaMethod( Problem* problem, double theta ) : TimeSolver( problem ), _theta( theta ) {}

void ThetaMethod::Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                       const blaze::DynamicMatrix<double>&,
                                                                       const blaze::DynamicMatrix<double>&,
                                                                       const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                           std::vector<blaze::DynamicMatrix<double>>& uNew,
                           std::vector<blaze::DynamicMatrix<double>>& u,
                           std::vector<blaze::DynamicMatrix<double>>& lambda ) {
    if( _theta == 0 ) {    // explicit Euler
        for( int j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
            uNew[j] = u[j] - ( _dt / _dx ) * ( fluxFunc( lambda[j - 1], lambda[j], lambda[j + 1], lambda[j + 2] ) -
                                               fluxFunc( lambda[j - 2], lambda[j - 1], lambda[j], lambda[j + 1] ) );
        }
    }
    else if( _theta == 0.5 ) {    // Crank-Nicolson
    }
    else if( _theta == 1 ) {    // implicit Euler
    }
}
