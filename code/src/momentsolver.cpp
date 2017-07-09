#include "momentsolver.h"

MomentSolver::MomentSolver(Problem* problem) : _problem(problem)
{
    _quad = new Hermite(_problem->GetNQuadPoints());
    _closure = new Closure(_problem);

    _dx = (_problem->_b-_problem->_a)/_problem->_nCells;
    for( int j = 0; j<_nCells+4; ++j){
        _x[j] = (j-2)*_dx;
    }
    _dt = _dx*cfl/uL;
    _nTimeSteps = _problem->_tEnd/_dt;
}


void MomentSolver::Solve(){
    double t = 0;
    vector u, uNew;
    std::vector<vector> uNew, u, lambda;
    _uL = 12;
    _uR = 3.0;

    // create solution fields
    uNew.resize(_problem->_nCells);
    u.resize(_problem->_nCells);
    lambda.resize(_problem->_nCells);
    for( int j = 0; j<_nCells+4; ++j){
        lambda[j].resize(_problem->_nMoments);
        u[j].resize(_problem->_nMoments);
        uNew[j].resize(_problem->_nMoments);
    }

    // Begin time loop
    while( t < tEnd ){
        // Modify moments into realizable direction
        for( int j = 2; j<_nCells+2; ++j ){
            u[j] = CalculateMoments(lambda[j]);
        }
        // Time Update Moments
        for( int j = 2; j<_nCells+2; ++j ){
            uNew[j] = u[j] - (dt/dx)*(numFlux(lambda[j],lambda[j+1])-numFlux(lambda[j-1],lambda[j]));
        }
        // Time Update dual variables
        for( int j = 2; j<_nCells+2; ++j ){
            lambda[j] = _closure->SolveClosure(uNew[j],lambda[j]);
        }
        t += dt;
    }
}

vector MomentSolver::numFlux(vector lambda1,vector lambda2){
    vector out = vector(_problem->_nMoments);
    for(int i = 0; i<_problem->_nMoments; ++i){
        out[i] = 0.0;
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

std::vector<vector> MomentSolver::SetupIC(){
    std::vector<vector> = vector(_problem->nMoments);
    vector xi = _quad->getNodes();
    vector w = _quad->getWeight();
    phiTile = _closure->GetPhiTilde();
    for( int i = 0; i<_problem->_nMoments; ++i ){
        for(int j = 0; j<_problem->nCells; ++j){
            out[j][i] = 0.0;
            for( int k = 0; k<_problem->_nQuadPoints; ++k ){
                out[j] += w(k)*IC(_x(j),xi(k))*phiTile(i,k);
            }
        }
    }
    return out;
}

double MomentSolver::IC(double x,double xi){
    double a = 1.0;
    double b = 2.0;
    if( x < a ){
        return uL;
    }else if(x>a && x<b){
        return _uR+(_uL-_uR)*(b-x)/(b-a);
    }else{
        return _uR;
    }
}
