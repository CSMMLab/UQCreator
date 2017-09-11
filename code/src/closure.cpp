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
            _phiTilde[k][i] = _phi[k][i]*(2.0*i+1.0);
        }
    }
    // calculate partial matrix for Hessian calculation
    _hPartial = std::vector<blaze::DynamicMatrix<double> >(_nQuadPoints, blaze::DynamicMatrix<double>(_nMoments,_nMoments,0.0));
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ){
        for( int i = 0; i < _nMoments; ++i ){
            for( int j = 0; j < _nMoments; ++j ){
                _hPartial[k] = _phiTilde[k]*blaze::trans(_phiTilde[k])*w[k];
            }
        }
    }
    double du = 0.0;
    _uMinus = 3.0-du;
    _uPlus = 12.0+du;
}

double Closure::UKinetic(double Lambda){
    double ePos = exp(Lambda);
    double eNeg = 1/ePos;
    if(Lambda > 0){
        return _uPlus/(eNeg+1.0)+_uMinus*eNeg/(1.0+eNeg);
    }else{
        return _uMinus/(ePos+1.0)+_uPlus*ePos/(1.0+ePos);
    }
}

blaze::DynamicVector<double> Closure::UKinetic(const blaze::DynamicVector<double>& Lambda){
    blaze::DynamicVector<double> y(Lambda.size(),0.0);
    blaze::DynamicVector<double> ePos = blaze::exp(Lambda);
    blaze::DynamicVector<double> eNeg = blaze::DynamicVector<double>(Lambda.size(), 1.0)/ePos ;
    for(unsigned int k = 0; k<Lambda.size();++k){
        if(Lambda[k] > 0){
            y[k] = _uPlus/(eNeg[k]+1.0)+_uMinus*eNeg[k]/(1.0+eNeg[k]);
        }
        else{
            y[k] = _uMinus/(ePos[k]+1.0)+_uPlus*ePos[k]/(1.0+ePos[k]);
        }
    }
    return y;
}

double Closure::DUKinetic(double Lambda){
    double ePos = exp(Lambda);
    double eNeg = 1/ePos;
    if(Lambda > 0){
        return eNeg*(_uPlus - _uMinus )/(exp(-2.0*Lambda)+2.0*eNeg+1.0);
    }
    else{
        return ePos*( _uPlus - _uMinus )/(1.0+2.0*ePos+exp(2.0*Lambda));
    }
}

blaze::DynamicVector<double> Closure::SolveClosure(const blaze::DynamicVector<double>& u, blaze::DynamicVector<double> lambda){
    int maxRefinements = 1000;
    blaze::DynamicVector<double> g = Gradient(lambda, u);
    if( blaze::length(g) < _problem->GetEpsilon() ){
        return lambda;
    }
    blaze::DynamicVector<double> dlambda = -g;
    blaze::DynamicMatrix<double> H = Hessian(lambda);
    blaze::posv( H, g, 'L');
    blaze::DynamicVector<double> lambdaNew = lambda-g;
    blaze::DynamicVector<double> dlambdaNew = Gradient(lambdaNew, u);
    for( int l=0; l<_problem->GetMaxIterations(); ++l ){
        double stepSize = 1.0;
        if(l!=0){
            g = Gradient(lambda, u);
            dlambda = -g;
            H = Hessian(lambda);
            blaze::posv( H, g, 'L');
            lambdaNew = lambda-stepSize*g;
            dlambdaNew = Gradient(lambdaNew, u);
        }
        int refinementCounter = 0;
        while( blaze::length(dlambda) < blaze::length(dlambdaNew)){
            stepSize *= 0.5;
            lambdaNew = lambda-stepSize*g;
            dlambdaNew = Gradient(lambdaNew, u);
            if( blaze::length(dlambdaNew) < _problem->GetEpsilon() ){
                return lambdaNew;
            }
            else if(++refinementCounter > maxRefinements){
                std::cerr<<"[ERROR]: Newton needed too many refinement steps!"<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        lambda = lambdaNew;
        if( blaze::length(dlambdaNew) < _problem->GetEpsilon() ){
            return lambdaNew;
        }
    }
    std::cerr<<"[ERROR]: Newton did not converge!"<<std::endl;
    exit(EXIT_FAILURE);
    return blaze::DynamicVector<double>(u.size(),-1.0);

}

double Closure::EvaluateLambda(const blaze::DynamicVector<double>& lambda, int k){
    double tmp = 0;
    for(int i = 0; i<_nMoments; ++i ){
        tmp += lambda[i]*_phiTilde[k][i];
    }
    return tmp;
}

blaze::DynamicVector<double> Closure::EvaluateLambda(const blaze::DynamicVector<double>& lambda){
    blaze::DynamicVector<double> tmp(_nQuadPoints,0.0);
    for(int k = 0; k<_nQuadPoints; ++k){
        for(int i = 0; i<_nMoments; ++i ){
            tmp[k] += lambda[i]*_phiTilde[k][i];
        }
    }
    return tmp;
}

double Closure::EvaluateLambda(const blaze::DynamicVector<double>& lambda, const blaze::DynamicVector<double>& xi, int k){
    double tmp = 0;
    for(int i = 0; i<_nMoments; ++i ){
        tmp += lambda[i]*_basis->Evaluate(i,xi[k])*(2.0*i+1.0);
    }
    return tmp;
}

blaze::DynamicVector<double> Closure::EvaluateLambda(const blaze::DynamicVector<double>& lambda, const blaze::DynamicVector<double>& xi){
    blaze::DynamicVector<double> tmp(xi.size(),0.0);
    for(unsigned int k = 0; k<xi.size(); ++k){
        for(int i = 0; i<_nMoments; ++i ){
            tmp[k] += lambda[i]*_basis->Evaluate(i,xi[k])*(2.0*i+1.0);
        }
    }
    return tmp;
}

blaze::DynamicVector<double> Closure::Gradient(const blaze::DynamicVector<double>& lambda, const blaze::DynamicVector<double>& u){
    blaze::DynamicVector<double> g(_nMoments, 0.0);
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ){
        g += UKinetic(blaze::inner( lambda, _phiTilde[k] ))*_phiTilde[k]*w[k];
    }
    return 0.5*g-u;
}

blaze::DynamicMatrix<double> Closure::Hessian(const blaze::DynamicVector<double>& lambda){
    blaze::DynamicMatrix<double> H(_nMoments,_nMoments, 0.0);
    for( int k = 0; k < _nQuadPoints; ++k ){
        H += _hPartial[k]*DUKinetic(blaze::inner( lambda, _phiTilde[k] ));
    }
    return 0.5*H;
}
