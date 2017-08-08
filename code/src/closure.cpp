#include "closure.h"

Closure::Closure(Problem *problem): _problem(problem), _nMoments(_problem->GetNMoments()), _nQuadPoints(_problem->GetNQuadPoints())
{
    // initialize classes
    _basis = new Legendre(_nMoments);
    _quad = new Legendre(_nQuadPoints);
    // calculate basis functions evaluated at the quadrature points
    _phi = std::vector<blaze::DynamicVector<double> >(_nQuadPoints, blaze::DynamicVector<double>(_nMoments,0.0));
    _phiTilde = std::vector<blaze::DynamicVector<double> >(_nQuadPoints, blaze::DynamicVector<double>(_nMoments,0.0));
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    for( int k = 0; k < _nQuadPoints; ++k ){
        for( int i = 0; i < _nMoments; ++i){
            _phi[k][i] = _basis->Evaluate(i,xi[k]);
            _phiTilde[k][i] = _phi[k][i]*(2.0*double(i)+1.0)/2.0;
        }
    }
    // calculate partial matrix for Hessian calculation
    _hPartial = std::vector<blaze::DynamicMatrix<double> >(_nQuadPoints, blaze::DynamicMatrix<double>(_nMoments,_nMoments,0.0));
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ){
        for( int i = 0; i < _nMoments; ++i ){
            for( int j = 0; j < _nMoments; ++j ){
                _hPartial[k] = _phi[k]*blaze::trans(_phi[k])*w[k];
            }
        }
    }
    _uMinus = 0.0;
    _uPlus = 20.0;
}

double Closure::UKinetic(double Lambda){
    if(Lambda > 0){
        return _uPlus/(exp(-Lambda)+1.0)+_uMinus*exp(-Lambda)/(1.0+exp(-Lambda));
    }else{
        return _uMinus/(exp(Lambda)+1.0)+_uPlus*exp(Lambda)/(1.0+exp(Lambda));
    }
}

double Closure::DUKinetic(double Lambda){
    if(Lambda > 0){
        //std::cout<<exp(-Lambda)*(_uPlus - _uMinus )/(exp(-2.0*Lambda)+2.0*exp(-Lambda)+1.0)<<std::endl;
        return exp(-Lambda)*(_uPlus - _uMinus )/(exp(-2.0*Lambda)+2.0*exp(-Lambda)+1.0);
    }else{
        //std::cout<<"Negative Lambda: "<<exp(Lambda)*( _uPlus - _uMinus )/(1.0+2.0*exp(Lambda)+exp(2.0*Lambda))<<std::endl;
        //std::cout<<"Term 1: "<< _uPlus - _uMinus <<std::endl;
        return exp(Lambda)*( _uPlus - _uMinus )/(1.0+2.0*exp(Lambda)+exp(2.0*Lambda));
    }
}

blaze::DynamicVector<double> Closure::SolveClosure(blaze::DynamicVector<double> u, blaze::DynamicVector<double> lambda){
    double stepSize = 0.5;
    blaze::DynamicVector<double> dlambda;

    blaze::DynamicVector<double> g = Gradient(lambda, u);
    for( int l=0; l<_problem->GetMaxIterations(); ++l ){
        blaze::DynamicMatrix<double> H = Hessian(lambda);
        dlambda = -g;
        std::cout<<H<<std::endl;
        blaze::posv( H, dlambda, 'L');
        blaze::DynamicVector<double> lambdaNew = lambda+stepSize*dlambda;
        blaze::DynamicVector<double> gNew = Gradient(lambdaNew, u);
        if( blaze::sqrLength(gNew) < _problem->GetEpsilon() ){
            return lambdaNew;
        }
        while( blaze::sqrLength(g) < blaze::sqrLength(gNew) ){
            stepSize = stepSize*0.5;
            lambdaNew = lambda+stepSize*dlambda;
            gNew = Gradient(lambdaNew, u);
        }
        lambda = lambdaNew;
        g = gNew;
        std::cout<<blaze::sqrLength(g)<<std::endl;
    }

/*
    blaze::DynamicVector<double> g = Gradient(lambda, u);
    std::cout<<"--------"<<std::endl<<blaze::sqrLength(g)<<std::endl;
    while( blaze::sqrLength(g) > _problem->GetEpsilon() ){
        blaze::DynamicMatrix<double> H = Hessian(lambda);
        dlambda = -g;
        //std::cout<<H<<std::endl;
        blaze::posv( H, dlambda, 'L');
        lambda = lambda+stepSize*dlambda;

        g = Gradient(lambda, u);
        std::cout<<blaze::sqrLength(g)<<std::endl;
    }*/
    exit(EXIT_FAILURE);
    return lambda;
}

double Closure::EvaluateLambda(blaze::DynamicVector<double> lambda, double xi){
    double tmp = 0;
    for(int i = 0; i<_nMoments-1; ++i ){
        tmp += lambda[i]*_basis->Evaluate(i,xi);
    }
    return tmp;
}

blaze::DynamicVector<double> Closure::Gradient(blaze::DynamicVector<double> lambda, blaze::DynamicVector<double> u){
    blaze::DynamicVector<double> g = blaze::DynamicVector<double>(_nMoments, 0.0);
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ){
        g += UKinetic(inner( lambda, _phi[k] ))*_phi[k]*w[k];
    }
   // std::cout<<g<<std::endl;
   // std::cout<<u<<std::endl;
    return g-u;
}

blaze::DynamicMatrix<double> Closure::Hessian(blaze::DynamicVector<double> lambda){
    blaze::DynamicMatrix<double> H = blaze::DynamicMatrix<double>(_nMoments,_nMoments, 0.0);
    for( int k = 0; k < _nQuadPoints; ++k ){
        H += _hPartial[k]*DUKinetic(inner( lambda, _phi[k] ));
        //std::cout<<DUKinetic(inner( lambda, _phi[k] ))<<std::endl;
        //std::cout<<DUKinetic(0)<<std::endl;
    }
    return H;
}
