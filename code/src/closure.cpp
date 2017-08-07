#include "closure.h"

Closure::Closure(Problem *problem): _problem(problem), _nMoments(_problem->GetNMoments()), _nQuadPoints(_problem->GetNQuadPoints())
{
    // initialize classes
    _basis = new Legendre(_nMoments);
    // calculate basis functions evaluated at the quadrature points
    _phi = std::vector<blaze::DynamicVector<double> >(_nQuadPoints, blaze::DynamicVector<double>(_nMoments,0.0));
    _phiTilde = std::vector<blaze::DynamicVector<double> >(_nQuadPoints, blaze::DynamicVector<double>(_nMoments,0.0));
    blaze::DynamicVector<double> xi = _basis->GetNodes();
    for( int k = 0; k < _nQuadPoints; ++k ){
        for( int i = 0; i < _nMoments-1; ++i){
            _phi[k][i] = _basis->Evaluate(i,xi[k]);
            _phiTilde[k][i] = _phi[k][i]*(2.0*i+1)/2;
        }
    }
    // calculate partial matrix for Hessian calculation
    _hPartial.resize(_nQuadPoints);
    blaze::DynamicVector<double> w = _basis->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ){
        _hPartial[k].resize(_nMoments,_nMoments);
        for( int i = 0; i < _nMoments-1; ++i ){
            for( int j = 0; j < _nMoments-1; ++j ){
                _hPartial[k] = _phi[k]*blaze::trans(_phi[k])*w[k];
            }
        }
    }
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
        return exp(-Lambda)*(_uPlus - _uMinus )/(exp(-2.0*Lambda)+2.0*exp(-Lambda)+1.0);
    }else{
        return exp(Lambda)*( _uPlus - _uMinus )/(1.0+2.0*exp(Lambda)+exp(2.0*Lambda));
    }
}

blaze::DynamicVector<double> Closure::SolveClosure(blaze::DynamicVector<double> u, blaze::DynamicVector<double> lambda){
    double stepSize = 0.5;
    blaze::DynamicVector<double> dlambda;
    for( int l=0; l<_problem->GetMaxIterations(); ++l ){
        blaze::DynamicVector<double> g = Gradient(lambda, u);
        while( blaze::sqrLength(g) < _problem->GetEpsilon() ){
            blaze::DynamicMatrix<double> H = Hessian(lambda);
            dlambda = -g;
            blaze::posv( H, dlambda, 'L');
            lambda = lambda+stepSize*dlambda;
        }
    }
    return lambda;
}

double Closure::EvaluateLambda(blaze::DynamicVector<double> lambda, double xi){
    double tmp = 0;
    for(int i = 0; i<_nMoments-1; ++i ){
        tmp = tmp+lambda[i]*_basis->Evaluate(i,xi);
    }
    return tmp;
}

blaze::DynamicVector<double> Closure::Gradient(blaze::DynamicVector<double> lambda, blaze::DynamicVector<double> u){
    blaze::DynamicVector<double> g = blaze::DynamicVector<double>(_nMoments);
    blaze::DynamicVector<double> w = _basis->GetWeights();
    for( int i = 0; i<_nMoments; ++i ){
        g[i] = 0;
    }
    for( int k = 0; k < _nQuadPoints; ++k ){
        g = g + UKinetic(inner( lambda, _phi[k] ))*_phi[k]*w[k];
    }
    return g;
}

blaze::DynamicMatrix<double> Closure::Hessian(blaze::DynamicVector<double> lambda){
    blaze::DynamicMatrix<double> H = blaze::DynamicMatrix<double>(_nMoments,_nMoments);
    for( int i = 0; i<_nMoments; ++i ){
        for( int j = 0; j<_nMoments; ++j ){
            H(i,j) = 0;
        }
    }
    for( int k = 0; k < _nQuadPoints; ++k ){
        H = H + _hPartial[k]*DUKinetic(inner( lambda, _phi[k] ));
    }
    return H;
}
