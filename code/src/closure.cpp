#include "closure.h"
typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;

Closure::Closure(Problem *problem): _problem(problem), _nMoments(_problem->GetNMoments()), _nQuadPoints(_problem->GetNQuadPoints())
{
    //
    // initialize classes
    _quadrature = new Quadrature(_problem);
    _basis = new BasisFunctions(_problem);

    // calculate basis functions evaluated at the quadrature points
    _phi.resize(_nQuadPoints);
    _phiTilde.resize(_nQuadPoints);
    for( int k = 0; k < _nQuadPoints; ++k ){
        _phi[k].resize(_nMoments);
        _phiTilde[k].resize(_nMoments);
        for( int i = 0; i < _nMoments-1; ++i){
            _phi[k](i) = _basis->Calculate(i,_quadrature->GetQuadPoint(k));
            _phiTilde[k](i) = _phi[k](i)*2.0/(2.0*(i-1)+1);
        }
    }

    // calculate partial matrix for Hessian calculation
    _hPartial.resize(_nQuadPoints);
    for( int k = 0; k < _nQuadPoints; ++k ){
        _hPartial[k].resize(_nMoments,_nMoments);
        for( int i = 0; i < _nMoments-1; ++i ){
            for( int j = 0; j < _nMoments-1; ++j ){
                _hPartial[k] = _phi[k]*transpose(_phi[k])*_quadrature->GetWeight(k);
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
    for( int l; l<_problem->GetMaxIterations(); ++l ){
        g = Gradient(lambda, u);
        while( blaze::sqrLength(g) < _problem->GetEpsilonNewton() ){
            vector dlambda = blaze::posv( Hessian(lambda), -g, 'L');
            lambda = lambda+dx*dlambda;
        }
    }
}

double Closure::EvaluateLambda(vector lambda,double xi){
    double tmp = 0;
    for(int i = 0; i<_nMoments-1; ++i ){
        tmp = tmp+lambda(i)*_basis->Calculate(i,xi);
    }
    return tmp;
}

vector Closure::EvaluateLambda(vector lambda,vector xi){
    vector tmp = vector(_nQuadPoints);
    for(int i = 0; i<_nMoments-1; ++i ){
        for( int k = 0; k<_nQuadPoints; ++k ){
            tmp[k] = tmp[k]+lambda(i)*_phi[k];
        }
    }
    return tmp;
}

vector Closure::Gradient(vector lambda, vector u){
    vector g = vector(_nMoments);
    for( int i = 0; i<_nMoments; ++i ){
        g(i) = 0;
    }
    for( int k = 0; k < _nQuadPoints; ++k ){
        g = g + UKinetic(inner( lambda, _phi[k] ))*_phi[k]*quadrature->GetWeight(k);
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
