#include "Rhs1D.h"
#include "MinMod1D.h"
#include "VanLeer1D.h"

Rhs1D::Rhs1D( Settings* settings, PhysicalProblem* physProb, double dx ) : Rhs( settings, physProb, dx ) {
    _riemannSolver = new HLLSolver1D( _physProb );

    if( _settings->GetLimiter() == VANLEER ) {
        _limiter = new VanLeer1D( _physProb );
    }
    else if( _settings->GetLimiter() == MINMOD ) {
        _limiter = new MinMod1D( _physProb );
    }
    else {
        _limiter = 0;
    }
    _t1_limitRightLeft = new Vector( _physProb->GetStateDim() );
    _t1_limitRightLeft->MakeZero();
    _limitRightLeft = new Vector( _physProb->GetStateDim() );

    _t1_limitLeftCurrent = new Vector( _physProb->GetStateDim() );
    _t1_limitLeftCurrent->MakeZero();
    _limitLeftCurrent = new Vector( _physProb->GetStateDim() );

    _t1_limitRightCurrent = new Vector( _physProb->GetStateDim() );
    _t1_limitRightCurrent->MakeZero();
    _limitRightCurrent = new Vector( _physProb->GetStateDim() );

    _t1_limitLeftRight = new Vector( _physProb->GetStateDim() );
    _t1_limitLeftRight->MakeZero();
    _limitLeftRight = new Vector( _physProb->GetStateDim() );
}

Rhs1D::~Rhs1D() {
    delete _riemannSolver;
    delete _limiter;
    delete _limitRightLeft;
    delete _limitLeftCurrent;
    delete _limitRightCurrent;
    delete _limitLeftRight;
}

void Rhs1D::SetupRhs( int index, Mesh* mesh, Vector* result ) {
    result->MakeZero();
    // numerical fluxes over right boundary, flux over left boundary is stored in result
    Vector* fluxRight = new Vector( _physProb->GetStateDim() );
    fluxRight->MakeZero();

    if( _settings->GetLimiter() != NOLIMITER ) {
        // calculate limited flux
        // if(index>0 && index<mesh->GetNumberOfCells()-1)
        //{
        _limiter->CalculateLimitedState1D( index - 1, mesh, _limitLeftCurrent, _limitRightLeft );    //_limitLeftCurrent works as dummy
        _limiter->CalculateLimitedState1D( index + 1, mesh, _limitLeftRight, _limitLeftCurrent );    //_limitLeftCurrent works as dummy
        _limiter->CalculateLimitedState1D( index, mesh, _limitLeftCurrent, _limitRightCurrent );
        /*}
        else if(index==0)
        {
            _limiter->CalculateLimitedState1D(index+1,mesh,_limitLeftRight,_limitLeftCurrent);//_limitLeftCurrent works as dummy
            _limiter->CalculateLimitedState1D(index,mesh,_limitLeftCurrent,_limitRightCurrent);
            //(*_limitRightLeft)=(*_limitLeftCurrent);
            _limiter->CalculateLimitedState1D(index-1,mesh,_limitLeftCurrent,_limitRightLeft);
        }
        else if(index==mesh->GetNumberOfCells()-1)
        {
            _limiter->CalculateLimitedState1D(index-1,mesh,_limitLeftCurrent,_limitRightLeft);//_limitLeftCurrent works as dummy
            _limiter->CalculateLimitedState1D(index,mesh,_limitLeftCurrent,_limitRightCurrent);
            //(*_limitLeftRight)=(*_limitRightCurrent);
            _limiter->CalculateLimitedState1D(index+1,mesh,_limitLeftRight,_limitLeftCurrent);
        }*/
    }
    else {
        ( *_limitRightLeft )    = ( *mesh->GetLeftCell( index )->GetState() );
        ( *_limitLeftCurrent )  = ( *mesh->GetCell( index )->GetState() );
        ( *_limitRightCurrent ) = ( *_limitLeftCurrent );
        ( *_limitLeftRight )    = ( *mesh->GetRightCell( index )->GetState() );
    }
    // calculation of fluxes with Riemann Solver
    _riemannSolver->NumericalFlux( _limitRightLeft, _limitLeftCurrent, result );
    _riemannSolver->NumericalFlux( _limitRightCurrent, _limitLeftRight, fluxRight );

    // result = fluxLeft-fluxRight
    result->Subtract( fluxRight );

    // result = (fluxLeft-fluxRight)/dx
    result->Multiply( 1 / _dx );

    delete fluxRight;
}

void Rhs1D::t1_SetupRhs( int index, Mesh* mesh, Mesh* t1_mesh, Vector* result, Vector* t1_result ) {
    t1_result->MakeZero();
    result->MakeZero();
    // numerical fluxes over right boundary, flux over left boundary is stored in result
    Vector* t1_fluxRight = new Vector( _physProb->GetStateDim() );
    Vector* fluxRight    = new Vector( _physProb->GetStateDim() );
    t1_fluxRight->MakeZero();
    fluxRight->MakeZero();

    if( _settings->GetLimiter() != NOLIMITER ) {
        // calculate limited flux
        // if(index>0 && index<mesh->GetNumberOfCells()-1)
        //{
        _limiter->t1_CalculateLimitedState1D(
            index - 1, mesh, t1_mesh, _limitLeftCurrent, _t1_limitLeftCurrent, _limitRightLeft, _t1_limitRightLeft );
        _limiter->CalculateLimitedState1D( index - 1, mesh, _limitLeftCurrent, _limitRightLeft );    //_limitLeftCurrent works as dummy

        _limiter->t1_CalculateLimitedState1D(
            index + 1, mesh, t1_mesh, _limitLeftRight, _t1_limitLeftRight, _limitLeftCurrent, _t1_limitLeftCurrent );
        _limiter->CalculateLimitedState1D( index + 1, mesh, _limitLeftRight, _limitLeftCurrent );    //_limitLeftCurrent works as dummy

        _limiter->t1_CalculateLimitedState1D(
            index, mesh, t1_mesh, _limitLeftCurrent, _t1_limitLeftCurrent, _limitRightCurrent, _t1_limitRightCurrent );
        _limiter->CalculateLimitedState1D( index, mesh, _limitLeftCurrent, _limitRightCurrent );
        /*}
        else if(index==0)
        {
            _limiter->t1_CalculateLimitedState1D(index+1,mesh,t1_mesh,_limitLeftRight,_t1_limitLeftRight,_limitLeftCurrent,_t1_limitLeftCurrent);
            _limiter->CalculateLimitedState1D(index+1,mesh,_limitLeftRight,_limitLeftCurrent);//_limitLeftCurrent works as dummy

            _limiter->t1_CalculateLimitedState1D(index,mesh,t1_mesh,_limitLeftCurrent,_t1_limitLeftCurrent,_limitRightCurrent,_t1_limitRightCurrent);
            _limiter->CalculateLimitedState1D(index,mesh,_limitLeftCurrent,_limitRightCurrent);

            (*_t1_limitRightLeft)=(*_t1_limitLeftCurrent);
            (*_limitRightLeft)=(*_limitLeftCurrent);
        }
        else if(index==mesh->GetNumberOfCells()-1)
        {

            _limiter->t1_CalculateLimitedState1D(index-1,mesh,t1_mesh,_limitLeftCurrent,_t1_limitLeftCurrent,_limitRightLeft,_t1_limitRightLeft);
            _limiter->CalculateLimitedState1D(index-1,mesh,_limitLeftCurrent,_limitRightLeft);//_limitLeftCurrent works as dummy

            _limiter->t1_CalculateLimitedState1D(index,mesh,t1_mesh,_limitLeftCurrent,_t1_limitLeftCurrent,_limitRightCurrent,_t1_limitRightCurrent);
            _limiter->CalculateLimitedState1D(index,mesh,_limitLeftCurrent,_limitRightCurrent);

            (*_t1_limitLeftRight)=(*_t1_limitRightCurrent);
            (*_limitLeftRight)=(*_limitRightCurrent);
        }*/
    }
    else {
        ( *_t1_limitRightLeft ) = ( *t1_mesh->GetLeftCell( index )->GetState() );
        ( *_limitRightLeft )    = ( *mesh->GetLeftCell( index )->GetState() );

        ( *_t1_limitLeftCurrent ) = ( *t1_mesh->GetCell( index )->GetState() );
        ( *_limitLeftCurrent )    = ( *mesh->GetCell( index )->GetState() );

        ( *_t1_limitRightCurrent ) = ( *_t1_limitLeftCurrent );
        ( *_limitRightCurrent )    = ( *_limitLeftCurrent );

        ( *_t1_limitLeftRight ) = ( *t1_mesh->GetRightCell( index )->GetState() );
        ( *_limitLeftRight )    = ( *mesh->GetRightCell( index )->GetState() );
    }
    // calculation of fluxes with Riemann Solver
    _riemannSolver->t1_NumericalFlux( _limitRightLeft, _t1_limitRightLeft, _limitLeftCurrent, _t1_limitLeftCurrent, result, t1_result );
    _riemannSolver->NumericalFlux( _limitRightLeft, _limitLeftCurrent, result );

    _riemannSolver->t1_NumericalFlux( _limitRightCurrent, _t1_limitRightCurrent, _limitLeftRight, _t1_limitLeftRight, fluxRight, t1_fluxRight );
    _riemannSolver->NumericalFlux( _limitRightCurrent, _limitLeftRight, fluxRight );

    // result = fluxLeft-fluxRight
    t1_result->Subtract( t1_fluxRight );
    result->Subtract( fluxRight );

    // result = (fluxLeft-fluxRight)/dx
    t1_result->Multiply( 1 / _dx );
    result->Multiply( 1 / _dx );

    delete t1_fluxRight;
    delete fluxRight;
}
