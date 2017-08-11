#include "hermite.h"

Hermite::Hermite(int degree):Polynomial(degree){
    _nodes.resize(degree);
    _weights.resize(degree);
    Compute();
}

void Hermite::Compute(){
    assert(_degree > 0);

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(_degree, _degree, 0.0);

    for(int i=0; i<_degree-1; ++i){
        CM(i+1,i) = std::sqrt((i+1)/2);
        CM(i,i+1) = std::sqrt((i+1)/2);
    }

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix(CM);

    for(int i=0; i<_degree; ++i){
        _nodes[i] = evSys.first[i];
        _weights[i] = std::pow(evSys.second(0,i),2) * std::sqrt(PI);
    }
    Sort();
}

double Hermite::Evaluate(int m, double x){
    return boost::math::hermite(m,x);
}

blaze::DynamicVector<double> Hermite::GetNodes(){
    return _nodes;
}

blaze::DynamicVector<double> Hermite::GetWeights(){
    return _weights;
}
