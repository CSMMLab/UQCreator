#include "burgerssolver.h"

BurgersSolver::BurgersSolver(int nCells, double tEnd, double cfl, double a, double b, double uL, double uR): _nCells(nCells), _a(a), _b(b)
{
    _u = new double[_nCells+4];
    _x = new double[_nCells+4];
    _dx = (_b-_a)/_nCells;
    for( int j = 0; j<_nCells+4; ++j){
        _x[j] = (j-2)*_dx;
        _u[j] = IC(_x[j],uL,uR);
    }
    _dt = _dx*cfl/uL;
    _nTimeSteps = tEnd/_dt;

}

double BurgersSolver::H(double u, double v, double w){
    return v-(_dt/_dx)*(g(v,w)-g(u,v));
}

double BurgersSolver::g(double u, double v){
    return f(u);
}

double BurgersSolver::f(double u){
    return 0.5*u*u;
}

void BurgersSolver::Solve(){
    double* uNew = new double[_nCells+4];
    for( int n = 0; n<_nTimeSteps; ++n){
        for( int j = 2; j<_nCells+2; ++j){
            uNew[j] = H(_u[j-1],_u[j],_u[j+1]);
        }

        for( int j = 2; j<_nCells+2; ++j){
            _u[j] = uNew[j];
        }
    }
}

double BurgersSolver::IC(double x, double uL, double uR){
    double a = 1.0;
    double b = 2.0;
    if( x < a ){
        return uL;
    }else if(x>a && x<b){
        return uR+(uL-uR)*(b-x)/(b-a);
    }else{
        return uR;
    }
}

void BurgersSolver::Print()const{
    std::ofstream out("outFile");
    for( int j = 2; j < _nCells+2; ++j){
        out<<_x[j]<<" "<<_u[j]<<std::endl;
    }
}

