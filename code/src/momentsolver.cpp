#include "momentsolver.h"

MomentSolver::MomentSolver(Problem* problem) : _problem(problem)
{
    _quad = new Legendre(_problem->GetNQuadPoints());
    _closure = new Closure(_problem);
    _mesh = _problem->GetMesh();

    _x = _mesh->GetGrid();
    _dx = _mesh->GetSpacing()[0]; // only equidistant!
    _nCells = _mesh->GetNumCells();
    _a = _x[0];
    _b = _x[_x.size()-1];
    std::cout<<_x.size()<<std::endl;
    std::cout<<"a = "<<_a<<",  b="<<_b<<std::endl;

    //std::cout<<_dx<<std::endl;
    //std::cout<<_x[1]-_x[0]<<std::endl;
    //exit(EXIT_FAILURE);

    _dt = _dx*_problem->GetCFL()/12.0;
    _nTimeSteps = _problem->GetTEnd()/_dt;
    _nMoments = _problem->GetNMoments();
    _tEnd = _problem->GetTEnd();
}

void MomentSolver::Solve(){
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    double t = 0;
    std::vector<blaze::DynamicVector<double>> uNew, u;
    _uL = 12;
    _uR = 3.0;
    // create solution fields
    uNew = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));
    u = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));
    _lambda = std::vector<blaze::DynamicVector<double> >(_nCells+4, blaze::DynamicVector<double>(_nMoments, 0.0));

    u = SetupIC();
    //std::cout<<"moment is "<<u[39]<<std::endl;
    _lambda[39] = _closure->SolveClosure(u[39],_lambda[39]);
    std::cout<<"lambda is "<<_lambda[39]<<std::endl;
    blaze::DynamicVector<double> gNew = _closure->Gradient( _lambda[39], u[39]);
    std::cout<<blaze::sqrLength(gNew)<<std::endl;
    //exit(EXIT_FAILURE);

    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    std::vector<blaze::DynamicVector<double>> phiTilde = _closure->GetPhiTilde();
    std::vector<blaze::DynamicVector<double>> phi = _closure->GetPhi();

    for(int j = 0; j<_nCells+4; ++j){
        _lambda[j] = _closure->SolveClosure(u[j],_lambda[j]);
    }

    //exit(EXIT_FAILURE);
    // Begin time loop
    while( t < _tEnd ){
        // Modify moments into realizable direction
        for( int j = 2; j<_nCells+2; ++j ){
            u[j] = CalculateMoments(_lambda[j]);
        }
        // Time Update Moments
        for( int j = 2; j<_nCells+2; ++j ){
            uNew[j] = u[j] - (_dt/_dx)*(numFlux(_lambda[j],_lambda[j+1])-numFlux(_lambda[j-1],_lambda[j]));
        }
        // Time Update dual variables
        for( int j = 2; j<_nCells+2; ++j ){
            std::cout<<"-----"<<std::endl;
            _lambda[j] = _closure->SolveClosure(uNew[j],_lambda[j]);
        }

        if(t == 0)
            std::cout << std::fixed << std::setprecision(8) << "t = " << t << std::flush;
        else
            std::cout << std::fixed << std::setprecision(8) << "\r" << "t = " << t << std::flush;
        t += _dt;
    }
    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    std::cout << "\nFinished!\nRuntime: " << std::setprecision(3) << std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()/1000.0 << "s" <<std::endl;
}

blaze::DynamicVector<double> MomentSolver::numFlux(blaze::DynamicVector<double> lambda1,blaze::DynamicVector<double> lambda2){
    blaze::DynamicVector<double> out(_nMoments,0.0);
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        out += w[k]*_problem->G(_closure->UKinetic(_closure->EvaluateLambda(lambda1,xi[k])), _closure->UKinetic(_closure->EvaluateLambda(lambda2,xi[k])))*_closure->GetPhiTilde(k);
    }
    return 0.5*out;
}

blaze::DynamicVector<double> MomentSolver::CalculateMoments(blaze::DynamicVector<double> lambda){
    blaze::DynamicVector<double> out(_nMoments,0.0);
    blaze::DynamicVector<double> xi(_quad->GetNodes());
    blaze::DynamicVector<double> w(_quad->GetWeights());
    for( int k = 0; k<_problem->GetNQuadPoints(); ++k){
        out += w[k]*_closure->UKinetic(_closure->EvaluateLambda(lambda,xi[k]))*_closure->GetPhiTilde(k);
    }
    return 0.5*out;
}

std::vector<blaze::DynamicVector<double>> MomentSolver::SetupIC(){
    std::vector<blaze::DynamicVector<double>> out(_nCells+4, blaze::DynamicVector<double>(_nMoments,0.0));
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w = _quad->GetWeights();
    std::vector<blaze::DynamicVector<double>> phiTilde = _closure->GetPhiTilde();
    for(int j = 0; j<_nCells+4; ++j){
        for( int i = 0; i<_nMoments; ++i ){
            for( int k = 0; k<_problem->GetNQuadPoints(); ++k ){
                out[j][i] += 0.5*w[k]*IC(_x[j],xi[k])*phiTilde[k][i];
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
    int cellIndex, nXi;
    try{
        auto file = cpptoml::parse_file(_problem->GetInputFile());
        auto plot = file->get_table("plot");

        cellIndex = plot->get_as<int>("cellIndex").value_or(-1);
        nXi = plot->get_as<int>("evalPoints").value_or(-1);
    }
    catch (const cpptoml::parse_exception& e){
        std::cerr << "Failed to parse " << _problem->GetInputFile() << ": " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    blaze::DynamicVector<double> xi(nXi,0.0);
    for( int k = 0; k<nXi; ++k ){
        xi[k] = -1 + k/(nXi-1.0)*2.0;
    }
    std::cout<<_x[cellIndex-1]<<std::endl;
    std::cout<<_lambda[cellIndex-1]<<std::endl;
    blaze::DynamicVector<double> res = _closure->UKinetic(_closure->EvaluateLambda(_lambda[cellIndex-1],xi));

    blaze::DynamicVector<double> test(2);
    blaze::DynamicVector<double> uTest(2);
    test[0] = -2.697521548829723;
    test[1] = 0.725439642287793;
    uTest[0] = 3.948033977800924;
    uTest[1] = 1.469034953656129;

    //std::cout<<"Gradient is "<<_closure->Gradient(test,uTest)<<std::endl;

    blaze::DynamicVector<double> exRes(nXi);
    for( int k = 0; k<nXi; ++k ){
        exRes[k] = _problem->ExactSolution(_tEnd,_x[cellIndex-1],xi[k]);
    }

    std::vector<double> xiStdVec(res.size(),0.0), resStdVec(res.size(),0.0), exResStdVec(res.size(),0.0);
    for(int i=0; i<nXi; i++){
        xiStdVec[i] = xi[i];
        resStdVec[i] = res[i];
        exResStdVec[i] = exRes[i];
    }

    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d(std::make_pair(xiStdVec, resStdVec)) << "with lines, " << gp.file1d(std::make_pair(xiStdVec, exResStdVec)) << "with lines\n";
}

void MomentSolver::PlotFixedXi(){
    double xi = 1.0;

    blaze::DynamicVector<double> res(_nCells);

    blaze::DynamicVector<double> exRes(_nCells);
    for( int k = 0; k<_nCells; ++k ){
        exRes[k] = _problem->ExactSolution(_tEnd,_x[k],xi);
        res[k] = _closure->UKinetic(_closure->EvaluateLambda(_lambda[k],xi));
    }


    std::vector<double> xStdVec(res.size(),0.0), resStdVec(res.size(),0.0), exResStdVec(res.size(),0.0);
    for(int i=0; i<_nCells; i++){
        xStdVec[i] = _x[i];
        resStdVec[i] = res[i];
        exResStdVec[i] = exRes[i];
    }

    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d(std::make_pair(xStdVec, resStdVec)) << "with lines, " << gp.file1d(std::make_pair(xStdVec, exResStdVec)) << "with lines\n";
}

void MomentSolver::Print(){
    int cellIndex = 55;
    int nXi = 20;
    blaze::DynamicVector<double> xi(nXi,0.0);
    for( int k = 0; k<nXi; ++k ){
        xi[k] = -1 + k/(nXi-1.0)*2.0;
    }
    std::cout<<_x[cellIndex]<<std::endl;
    std::cout<<_closure->UKinetic(_closure->EvaluateLambda(_lambda[cellIndex],xi))<<std::endl;
}
