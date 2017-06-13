#include "hermite.h"

void Hermite::computeNodes(int degree){

    //construct companion matrix
    blaze::DynamicMatrix<double> CM(degree, degree);

    for(double i=0; i<degree-1; ++i){
        CM(i+1,i) = std::sqrt((i+1)/2);
        CM(i,i+1) = std::sqrt((i+1)/2);
    }

    auto evSys = MathTools::computeEigenValTriDiagMatrix(CM);

    for(int i=0; i<degree; ++i){
        _nodes(i) = evSys.first[i]*std::sqrt(2);
        _weights(i) = std::pow(evSys.second[0][i],2);
    }
}
