#include "sspmultistep.h"

SSPMultiStep::SSPMultiStep( Problem* problem, Closure* closure ) : TimeSolver( problem ) {
    _counter = 0;
    _heun    = new Heun( problem, closure );
}

SSPMultiStep::~SSPMultiStep() { delete _heun; }

void SSPMultiStep::Advance( std::function<blaze::DynamicMatrix<double>( const blaze::DynamicMatrix<double>&,
                                                                        const blaze::DynamicMatrix<double>&,
                                                                        const blaze::DynamicMatrix<double>&,
                                                                        const blaze::DynamicMatrix<double>& )> const& fluxFunc,
                            std::vector<blaze::DynamicMatrix<double>>& uNew,
                            std::vector<blaze::DynamicMatrix<double>>& u,
                            std::vector<blaze::DynamicMatrix<double>>& lambda ) {
    _counter++;
    if( _counter < 4 ) {
        _heun->Advance( fluxFunc, uNew, u, lambda );
    }
    else {
        for( int j = 3; j < _problem->GetMesh()->GetNumCells() + 1; ++j ) {
            uNew[j] = ( 8.0 / 9.0 ) * u[j] + ( 1.0 / 9.0 ) * _u3Step[j] -
                      ( 4.0 / 3.0 ) * ( _dt / _dx ) *
                          ( fluxFunc( lambda[j - 1], lambda[j], lambda[j + 1], lambda[j + 2] ) -
                            fluxFunc( lambda[j - 2], lambda[j - 1], lambda[j], lambda[j + 1] ) );
        }
    }
    _u3Step = _u2Step;
    _u2Step = _u1Step;
    _u1Step = u;
}
