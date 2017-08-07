#include "legendre.h"

Legendre::Legendre(int degree):Polynomial(degree){
    _nodes = blaze::DynamicVector<double>(degree,0.0);
    _weights = blaze::DynamicVector<double>(degree,0.0);
    ComputeNodes(degree);
}

void Legendre::ComputeNodes(int degree){

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(degree, degree, 0.0);

    for(double i=0; i<degree-1; ++i){
        CM(i+1,i) = std::sqrt(1/(4-1/std::pow(i+1,2)));
        CM(i,i+1) = std::sqrt(1/(4-1/std::pow(i+1,2)));
    }

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix(CM);

    for(int i=0; i<degree; ++i){
        if(std::fabs(evSys.first[i])<1e-15)
            _nodes[i] = 0;
        else
            _nodes[i] = evSys.first[i];
        _weights[i] = 2*std::pow(evSys.second(0,i),2);
    }
    Sort();
}

double Legendre::Evaluate(int m,double x){
    return boost::math::legendre_p(m, x);
}

blaze::DynamicVector<double> Legendre::GetNodes(){
    return _nodes;
}

blaze::DynamicVector<double> Legendre::GetWeights(){
    return _weights;
}
