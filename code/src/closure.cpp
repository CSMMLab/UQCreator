#include "closure.h"

Closure::Closure(Problem problem): _problem(problem)
{
    // initialize classes
    _newton = new Newton(_problem);
    _quadrature = new Quadrature(_problem);
    _basis = new BasisFunctions(_problem);

    // initialize basis functions evaluated at the quadrature points
    _phi = new blaze::DynamicMatrix<float,blaze::rowMajor>(_problem->GetNMoments(),_problem->GetNQuadPoints());

    for( int i = 0; i < _problem->GetNMoments()-1; ++i){
        for( int k = 0; k < _problem->GetNQuadPoints()-1; ++k ){
            (*_phi)(i,k) = _basis->Calculate(_quadrature->GetQuadpoint(k));
        }
    }

    // calculate partial matrix for Hessian calculation
    _hPartial = new vector< blaze::DynamicMatrix<float>* >(_problem->GetNQuadPoints());
    for( int k = 0; k < _problem->GetNQuadPoints(); ++k ){
        (*_hPartial)[k] = new blaze::DynamicMatrix<float>(_problem->GetNMoments(),_problem->GetNMoments());
    }


    for( int k = 0; k < _problem->GetNQuadPoints()-1; ++k ){
        for( int i = 0; i < _problem->GetNMoments()-1; ++i ){
            for( int j = 0; j < _problem->GetNMoments()-1; ++j ){
                (*_hPartial->at(k))(i,j) = (*_phi)(i,k)*(*_phi)(j,k)*_quadrature->GetWeight(k);
            }
        }
    }

}

blaze::DynamicVector<double,blaze::rowVector> Closure::SolveClosure(blaze::DynamicVector<double,blaze::rowVector> u){
    Math::Newton();
}
