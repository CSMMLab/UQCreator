#include "momentsolver.h"

MomentSolver::MomentSolver(Problem* problem) : _problem(problem)
{
    _quad = new Legendre(_problem->GetNQuadPoints());
    _closure = new Closure(_problem);
    _mesh = _problem->GetMesh();

    _x = _mesh->GetGrid();
    _dx = _mesh->GetSpacing()[0]; // only aequidistant!
    _nCells = _mesh->GetNumCells();
    _a = _x[0];
    _b = _x[_nCells-1];

    _dt = _dx*_problem->GetCFL()/12.0;
    _nTimeSteps = _problem->GetTEnd()/_dt;
    _nMoments = _problem->GetNMoments();
    _tEnd = _problem->GetTEnd();
}

void MomentSolver::Solve(){
    double t = 0;
    std::vector<blaze::DynamicVector<double>> uNew, u;
    _uL = 12;
    _uR = 3.0;
    // create solution fields
    uNew = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));
    u = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));
    _lambda = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));

    u = SetupIC();

    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    std::vector<blaze::DynamicVector<double>> phiTilde = _closure->GetPhi();

    //std::cout << xi << std::endl;
    //std::cout << w << std::endl;

    double out = 0;
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k ){
        //out += w[k]*phiTilde[k][0]*phiTilde[k][2]*IC(_x[0],xi[k]);
        //std::cout<<phiTilde[k][4]<<std::endl;
        out += w[k]*phiTilde[k][4];//*IC(_x[8],xi[k]);
    }
    std::cout<<"out = "<<out<<std::endl;

    // DEBUG IC
    for(int j = 0; j<_nCells+4; ++j){
        std::cout << u[j] << std::endl;
    }

    for(int j = 0; j<_nCells+4; ++j){
        std::cout<<"Cell "<<j<<std::endl;
        _lambda[j] = _closure->SolveClosure(u[j],_lambda[j]);
        //std::cout << _lambda[j] << std::endl;
    }

    //exit(EXIT_FAILURE);
    // Begin time loop
    while( t < _tEnd ){
        // Modify moments into realizable direction
        for( int j = 2; j<_nCells+2; ++j ){
            std::cout<<"Cell "<<j<<std::endl;
            u[j] = CalculateMoments(_lambda[j]);
        }
        std::cout<<"Moments calculated for inexact dual states."<<std::endl;
        // Time Update Moments
        for( int j = 2; j<_nCells+2; ++j ){
            uNew[j] = u[j] - (_dt/_dx)*(numFlux(_lambda[j],_lambda[j+1])-numFlux(_lambda[j-1],_lambda[j]));
        }
        std::cout<<"Time Update calculated."<<std::endl;
        // Time Update dual variables
        for( int j = 2; j<_nCells+2; ++j ){
            _lambda[j] = _closure->SolveClosure(uNew[j],_lambda[j]);
        }
        std::cout<<"Dual states updated."<<std::endl;
        if(t == 0)
            std::cout << std::fixed << std::setprecision(8) << "t = " << t << std::flush;
        else
            std::cout << std::fixed << std::setprecision(8) << "\r" << "t = " << t << std::flush;
        t += _dt;
    }
    std::cout << std::endl;
}

blaze::DynamicVector<double> MomentSolver::numFlux(blaze::DynamicVector<double> lambda1,blaze::DynamicVector<double> lambda2){
    blaze::DynamicVector<double> out(_nMoments,0.0);
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        out += w[k]*_problem->G(_closure->UKinetic(_closure->EvaluateLambda(lambda1,xi[k])), _closure->UKinetic(_closure->EvaluateLambda(lambda2,xi[k])))*_closure->GetPhiTilde(k);
    }
    return out;
}

blaze::DynamicVector<double> MomentSolver::CalculateMoments(blaze::DynamicVector<double> lambda){
    blaze::DynamicVector<double> out(_nMoments,0.0);
    blaze::DynamicVector<double> xi(_quad->GetNodes());
    blaze::DynamicVector<double> w(_quad->GetWeights());
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        std::cout<<k<<std::endl;
        out += w[k]*_closure->UKinetic(_closure->EvaluateLambda(lambda,xi[k]))*_closure->GetPhiTilde(k);
    }
    return out;
}

std::vector<blaze::DynamicVector<double>> MomentSolver::SetupIC(){
    std::vector<blaze::DynamicVector<double>> out(_nCells+4, blaze::DynamicVector<double>(_nMoments,0.0));
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    std::vector<blaze::DynamicVector<double>> phiTilde = _closure->GetPhiTilde();
    for(int j = 0; j<_nCells+4; ++j){
        for( int i = 0; i<_nMoments; ++i ){
            for( int k = 0; k<_problem->GetNQuadPoints(); ++k ){
                out[j][i] += w[k]*IC(_x[j],xi[k])*phiTilde[k][i];
            }
        }
    }
    return out;
}

double MomentSolver::IC(double x,double xi){
    double a = 0.5;
    double b = 1.5;
    if( x < a+0.2*xi ){
        return _uL;
    }
    else if( x < b+0.2*xi ){
        return _uL+(_uR-_uL)*(a+0.2*xi-x)/(a-b);
    }
    else{
        return _uR;
    }
}

void MomentSolver::Plot(){

}

void MomentSolver::Print(){

}
