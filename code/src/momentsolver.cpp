#include "momentsolver.h"

MomentSolver::MomentSolver(Problem* problem) : _problem(problem)
{
    _quad = new Hermite(_problem->GetNQuadPoints());
    _closure = new Closure(_problem);

    _origSolver = new Burgers();

    _a = 0.0;
    _b = 3.0;
    _dx = (_b-_a)/_problem->_nCells;
    for( int j = 0; j<_nCells+4; ++j){
        _x[j] = (j-2)*_dx;
    }
    _dt = _dx*_problem->_cfl/12.0;
    _nTimeSteps = _origSolver->_tEnd/_dt;
    _nCells = _origSolver->_nCells;
    _nMoments = _problem->_nMoments;
    _tEnd = _origSolver->_tEnd;
}

void MomentSolver::Solve(){
    double t = 0;
    std::vector<vector> uNew, u, lambda;
    _uL = 12;
    _uR = 3.0;

    // create solution fields
    uNew.resize(_nCells);
    u.resize(_nCells);
    lambda.resize(_nCells);
    for( int j = 0; j<_nCells+4; ++j){
        lambda[j].resize(_nMoments);
        u[j].resize(_nMoments);
        uNew[j].resize(_nMoments);
    }

    // Begin time loop
    while( t < _tEnd ){
        // Modify moments into realizable direction
        for( int j = 2; j<_nCells+2; ++j ){
            u[j] = CalculateMoments(lambda[j]);
        }
        // Time Update Moments
        for( int j = 2; j<_nCells+2; ++j ){
            uNew[j] = u[j] - (_dt/_dx)*(numFlux(lambda[j],lambda[j+1])-numFlux(lambda[j-1],lambda[j]));
        }
        // Time Update dual variables
        for( int j = 2; j<_nCells+2; ++j ){
            lambda[j] = _closure->SolveClosure(uNew[j],lambda[j]);
        }
        t += _dt;
    }
}

vector MomentSolver::numFlux(vector lambda1,vector lambda2){
    vector out = vector(_nMoments);
    for(int i = 0; i<_nMoments; ++i){
        out[i] = 0.0;
    }
    vector xi = _quad->GetNodes();
    vector w = _quad->GetWeights();
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        out = out + w[k]*_origSolver->G(_closure->UKinetic(_closure->EvaluateLambda(lambda1,xi[k])), _closure->UKinetic(_closure->EvaluateLambda(lambda2,xi[k])))*_closure->GetPhiTilde(k);
    }
    return out;
}

vector MomentSolver::CalculateMoments(vector lambda){
    vector out = vector(_nMoments);
    for(int j = 0; j<_nCells; ++j){
        out[j] = 0.0;
    }
    vector xi = _quad->GetNodes();
    vector w = _quad->GetWeights();
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        out = out + w[k]*_closure->UKinetic(_closure->EvaluateLambda(lambda,xi[k]))*_closure->GetPhiTilde(k);
    }
    return out;
}

std::vector<vector> MomentSolver::SetupIC(){
    std::vector<vector> out;
    vector xi = _quad->GetNodes();
    vector w = _quad->GetWeights();
    std::vector<vector> phiTile = _closure->GetPhiTilde();
    for( int i = 0; i<_nMoments; ++i ){
        for(int j = 0; j<_nCells; ++j){
            out[j][i] = 0.0;
            for( int k = 0; k<_problem->GetNQuadPoints(); ++k ){
                out[j][i] = out[j][i] + w[k]*IC(_x[j],xi[k])*phiTile[i][k];
            }
        }
    }
    return out;
}

double MomentSolver::IC(double x,double xi){
    double a = 1.0;
    double b = 2.0;
    if( x < a ){
        return _uL;
    }else if(x>a && x<b){
        return _uR+(_uL-_uR)*(b-x)/(b-a);
    }else{
        return _uR;
    }
}
