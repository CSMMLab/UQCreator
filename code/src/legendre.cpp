#include "legendre.h"

#include <boost/math/special_functions/legendre.hpp>

Legendre::Legendre(int degree):Polynomial(degree){
    _nodes = vector(degree);
    _weights = vector(degree);
    computeNodes(degree);
}

void Legendre::computeNodes(int degree){

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(degree, degree);

    for(double i=0; i<degree-1; ++i){
        CM(i+1,i) = std::sqrt(1/(4-1/std::pow(i+1,2)));
        CM(i,i+1) = std::sqrt(1/(4-1/std::pow(i+1,2)));
    }

    auto evSys = MathTools::computeEigenValTriDiagMatrix(CM);

    for(int i=0; i<degree; ++i){
        _nodes[i] = evSys.first[i];
        _weights[i] = 2*std::pow(evSys.second(0,i),2);
    }
}

double Legendre::evaluate(int m,double x){
    return boost::math::legendre_p(m, x);
}

vector Legendre::getNodes(){
    return _nodes;
}

vector Legendre::getWeights(){
    return _weights;
}
