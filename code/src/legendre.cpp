#include "legendre.h"

#include <boost/math/special_functions/legendre.hpp>

Legendre::Legendre(int degree):Polynomial(degree){
    _nodes = vector(degree);
    _weights = vector(degree);
    ComputeNodes(degree);
}

void Legendre::ComputeNodes(int degree){

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(degree, degree);

    for(double i=0; i<degree-1; ++i){
        CM(i+1,i) = std::sqrt(1/(4-1/std::pow(i+1,2)));
        CM(i,i+1) = std::sqrt(1/(4-1/std::pow(i+1,2)));
    }

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix(CM);

    for(int i=0; i<degree; ++i){
        _nodes[i] = evSys.first[i];
        _weights[i] = 2*std::pow(evSys.second(0,i),2);
    }
}

double Legendre::Evaluate(int m,double x){
    return boost::math::legendre_p(m, x);
}

vector Legendre::GetNodes(){
    return _nodes;
}

vector Legendre::GetWeights(){
    return _weights;
}
