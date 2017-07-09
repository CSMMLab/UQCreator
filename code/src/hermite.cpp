#include "hermite.h"

Hermite::Hermite(int degree):Polynomial(degree){
    _nodes = vector(degree);
    _weights = vector(degree);
    ComputeNodes(degree);
}

void Hermite::ComputeNodes(int degree){

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(degree, degree);

    for(double i=0; i<degree-1; ++i){
        CM(i+1,i) = std::sqrt((i+1)/2);
        CM(i,i+1) = std::sqrt((i+1)/2);
    }

    auto evSys = MathTools::ComputeEigenValTriDiagMatrix(CM);

    for(int i=0; i<degree; ++i){
        _nodes[i] = evSys.first[i];
        _weights[i] = std::pow(evSys.second(0,i),2) * std::sqrt(PI);
    }
}

double Hermite::Evaluate(int m, double x){
    return -1.0;
}

vector Hermite::GetNodes(){
    return _nodes;
}

vector Hermite::GetWeights(){
    return _weights;
}
