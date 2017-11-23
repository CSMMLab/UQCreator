#include "thetamethod.h"

ThetaMethod::ThetaMethod(double theta) : _theta(theta){

}

double ThetaMethod::Solve(const blaze::DynamicVector<double>& u, const blaze::DynamicVector<double>& flux1, const blaze::DynamicVector<double>& flux2){
    blaze::DynamicVector<double> uNew(u.size(), 0.0);
    if(_theta == 0){            // explicit Euler
        //uNew[j] = u[j] - (_dt/_dx)*(numFlux(lambda[j-1], lambda[j], lambda[j+1], lambda[j+2]) - numFlux(lambda[j-2],_lambda[j-1],_lambda[j],_lambda[j+1]));
    }
    else if(_theta == 0.5){     // Crank-Nicolson

    }
    else if(_theta == 1){       // implicit Euler

    }
}
