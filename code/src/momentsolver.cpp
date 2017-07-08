#include "momentsolver.h"

MomentSolver::MomentSolver(Ploblem* problem): _problem(problem)
{
    _quad = new Hermite(_problem);
    _closure = new Closure(_problem);
    _origSolver = new BurgersSolver(nCells, tEnd, cfl, a, b, uL, uR);
}


void MomentSolver::Solve(){
    double t = 0;
    vector u, uNew;
    std::vector<vector> uNew, u, lambda;
    double uL = 12;
    double uR = 3.0;

    _dx = (_b-_a)/_problem->_nCells;
    uNew.resize(_problem->_nCells);
    u.resize(_problem->_nCells);
    lambda.resize(_problem->_nCells);
    for( int j = 0; j<_nCells+4; ++j){
        _x[j] = (j-2)*_dx;
        _u[j] = IC(_x[j],uL,uR);
    }
    _dt = _dx*cfl/uL;
    _nTimeSteps = _problem->_tEnd/_dt;

    // Begin time loop
    while( t < tEnd ){
        // Modify moments into realizable direction
        for( int j = 0; j<nCells; ++j ){
            u[j] = CalculateMoments(lambda[j]);
        }
        // Time Update Moments
        for( int j = 0; j<nCells; ++j ){
            uNew[j] = u[j] - (dt/dx)*(numFlux(lambda[j],lambda[j+1])-numFlux(lambda[j-1],lambda[j]));
        }
        // Time Update dual variables
        for( int j = 0; j<nCells; ++j){
            lambda[j] = _closure->SolveClosure(uNew[j],lambda[j]);
        }
        t += dt;
    }
}

vector MomentSolver::umFlux(vector lambda1,vector lambda2){
    vector out = vector(_problem->nMoments);
    for(int j = 0; j<nCells; ++j){
        out[j] = 0.0;
    }
    vector xi = _quad->getNodes();
    vector w = _quad->getWeight();
    for( int k = 0; k<_problem->nQuadPoints; ++k){
        out = out + w(k)*_origSolver->g(_closure->UKinetic(_closure->EvaluateLambda(lambda1,xi(k))), _closure->UKinetic(_closure->EvaluateLambda(lambda2,xi(k))))*_closure->GetPhiTilde(k);
    }
    return out;
}

vector MomentSolver::CalculateMoments(vector lambda){
    vector out = vector(_problem->nMoments);
    for(int j = 0; j<nCells; ++j){
        out[j] = 0.0;
    }
    vector xi = _quad->getNodes();
    vector w = _quad->getWeight();
    for( int k = 0; k<_problem->nQuadPoints; ++k){
        out = out + w(k)*_closure->UKinetic(_closure->EvaluateLambda(lambda1,xi(k)))*_closure->GetPhiTilde(k);
    }
    return out;
}
