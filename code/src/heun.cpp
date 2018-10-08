#include "heun.h"

Heun::Heun( Problem* _problem, Closure* closure ) : TimeSolver( _problem ), _closure( closure ) {}

Heun::~Heun() {}

void Heun::Advance( std::function<Matrix( const Matrix&,
                                                                const Matrix&,
                                                                const Matrix&,
                                                                const Matrix& )> const& fluxFunc,
                    std::vector<Matrix>& uNew,
                    std::vector<Matrix>& u,
                    std::vector<Matrix>& lambda ) {

    for( unsigned j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
        uNew[j] = u[j] - 0.5 * ( _dt / _dx ) *
                             ( fluxFunc( lambda[j - 1], lambda[j], lambda[j + 1], lambda[j + 2] ) -
                               fluxFunc( lambda[j - 2], lambda[j - 1], lambda[j], lambda[j + 1] ) );
    }

    for( unsigned j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
        lambda[j] = _closure->SolveClosure( uNew[j], lambda[j] );
    }

    for( unsigned j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
        uNew[j] = uNew[j] - 0.5 * ( _dt / _dx ) *
                                ( fluxFunc( lambda[j - 1], lambda[j], lambda[j + 1], lambda[j + 2] ) -
                                  fluxFunc( lambda[j - 2], lambda[j - 1], lambda[j], lambda[j + 1] ) );
    }
}
