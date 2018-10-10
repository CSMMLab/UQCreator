#include "Heun.h"
#include "Mesh1D.h"

Heun::Heun( Settings* settings, Mesh* mesh, PhysicalProblem* physicalProblem, Rhs* rhs )
    : TimeIntegrator( settings, mesh, rhs ), _physicalProblem( physicalProblem ) {
    _intermediates = new Mesh1D( _mesh->GetXLeft(), _mesh->GetXRight(), _mesh->GetDimX(), _physicalProblem, _settings );
}

Heun::~Heun() { delete _intermediates; }

void Heun::Integrate() {

    // calculate intermediates for all cells
    for( int i = 0; i < _intermediates->GetNumberOfCells(); ++i ) {
        _physicalProblem->CFL( _mesh->GetCell( i )->GetState(), _settings->GetDt(), _mesh->GetDx() );
        CalcIntermediate( i, _intermediates->GetCell( i )->GetState() );
    }

    // calculate next time step for all cells
    for( int i = 0; i < _intermediates->GetNumberOfCells(); ++i ) {
        UpdateWithIntermediate( i );
    }
}

void Heun::CalcIntermediate( int i, Vector* intermediate ) {
    // setting up right hand side for Vector current
    _rhs->SetupRhs( i, _mesh, intermediate );

    // intermediate = _dt*f(cell_current)
    intermediate->Multiply( _settings->GetDt() );

    // intermediate = cell_current+_dt*f(cell_current)
    intermediate->Add( _mesh->GetCell( i )->GetState() );
}

void Heun::UpdateWithIntermediate( int index ) {
    Vector update( _physicalProblem->GetStateDim() );

    // rhs with intermediate values
    _rhs->SetupRhs( index, _intermediates, &update );
    // update = _dt*f(x_intermediate)
    update.Multiply( _settings->GetDt() );

    // update = x_intermediate+_dt*f(x_intermediate)
    update.Add( _intermediates->GetCell( index )->GetState() );

    // update = x+x_intermediate+_dt*f(x_intermediate)
    update.Add( _mesh->GetCell( index )->GetState() );

    // update = 0.5*x+0.5*(x_intermediate+_dt*f(x_intermediate))
    update.Multiply( 0.5 );

    // store time update in cell
    _mesh->GetCell( index )->SetState( &update );
}
