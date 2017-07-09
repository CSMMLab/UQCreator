#include "closure.h"
typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;

Closure::Closure(Problem *problem): _problem(problem), _nMoments(_problem->GetNMoments()), _nQuadPoints(_problem->GetNQuadPoints())
{
    //
    // initialize classes
    _quadrature = new Quadrature(_problem);
    _basis = new Legendre(_nMoments);

    // calculate basis functions evaluated at the quadrature points
    _phi.resize(_nQuadPoints);
    _phiTilde.resize(_nQuadPoints);
    vector xi = _quadrature->GetNodes();
    for( int k = 0; k < _nQuadPoints; ++k ){
        _phi[k].resize(_nMoments);
        _phiTilde[k].resize(_nMoments);
        for( int i = 0; i < _nMoments-1; ++i){
            _phi[k][i] = _basis->Evaluate(i,xi[k]);
            _phiTilde[k][i] = _phi[k][i]*2.0/(2.0*(i-1)+1);
        }
    }

    // calculate partial matrix for Hessian calculation
    _hPartial.resize(_nQuadPoints);
    vector w = _quadrature->GetWeights();
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

vector Closure::SolveClosure(vector u, vector lambda){
    double stepSize = 0.5;
    vector dlambda;
    for( int l; l<_problem->GetMaxIterations(); ++l ){
        vector g = Gradient(lambda, u);
        while( blaze::sqrLength(g) < _problem->GetEpsilon() ){
            matrix H = Hessian(lambda);
            vector dlambda = -g;
            blaze::posv( H, dlambda, 'L');
            lambda = lambda+stepSize*dlambda;
        }
    }
}

double Closure::EvaluateLambda(vector lambda,double xi){
    double tmp = 0;
    for(int i = 0; i<_nMoments-1; ++i ){
        tmp = tmp+lambda[i]*_basis->Evaluate(i,xi);
    }
    return tmp;
}

vector Closure::Gradient(vector lambda, vector u){
    vector g = vector(_nMoments);
    vector w = _quadrature->GetWeights();
    for( int i = 0; i<_nMoments; ++i ){
        g[i] = 0;
    }
    for( int k = 0; k < _nQuadPoints; ++k ){
        g = g + UKinetic(inner( lambda, _phi[k] ))*_phi[k]*w[k];
    }
    return g;
}

matrix Closure::Hessian(vector lambda){
    matrix H = matrix(_nMoments,_nMoments);
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
