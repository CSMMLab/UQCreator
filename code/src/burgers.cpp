#include "burgers.h"

Burgers::Burgers(std::string inputFile) : Problem(inputFile){
    _x = _mesh->GetGrid();
    _dx = _mesh->GetSpacing().at(0);
    _dt = _dx*_CFL/12.0;
    _nCells = _mesh->GetNumCells();
    _nTimeSteps = _tEnd/_dt;
    _u = blaze::DynamicVector<double>(_nCells+4, 0.0);
}

double Burgers::H(double u, double v, double w){
    return v-(_dt/_dx)*(G(v,w)-G(u,v));
}

double Burgers::G(double u, double v){
    return F(u);
}

double Burgers::F(double u){
    return 0.5*u*u;
}

void Burgers::Solve(){
    double* uNew = new double[_nCells+4];
    for( int n = 0; n<_nTimeSteps; ++n){
        for( int j = 2; j<_nCells+2; ++j){
            //_dt = _timeDiscretization->getDt();
            uNew[j] = H(_u[j-1],_u[j],_u[j+1]);
        }

        for( int j = 2; j<_nCells+2; ++j){
            _u[j] = uNew[j];
        }
    }
}

double Burgers::IC(double x, double uL, double uR){
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

void Burgers::Print() const{
    std::ofstream out("outFile");
    for( int j = 2; j < _nCells+2; ++j){
        out<<_x[j]<<" "<<_u[j]<<std::endl;
    }
}

void Burgers::Plot(blaze::DynamicVector<double>& x, blaze::DynamicVector<double>& u) const{
    std::vector<double> x1, u1;
    for(int i=0; i<_nCells+4; i++){
        x1.push_back(x[i]);
        u1.push_back(u[i]);
    }
    Gnuplot gp;
    gp << "plot '-' with lines notitle\n";
    gp << "set xlabel 'Space'\n";
    gp << "set ylabel 'Velocity\n";
    gp.send1d(std::make_pair(x1, u1));
}

void Burgers::WriteToFile(std::string filename, int filetype) const{

}

double Burgers::ExactSolution(double t, double x, double xi){
    double uL = 12.0;
    double uR = 3.0;
    double sigma = 0.2;
    double x0 = 0.5;
    double x1 = 1.5;
    double x0Bar = x0+sigma*xi+uL*t;
    double x1Bar = x1+sigma*xi+uR*t;
    //double x0Bar = x0+sigma*xi+uL*(1.0/9.0)+0.5*(uL+uR)*(t-(1.0/9.0));
    //double x1Bar = x1+sigma*xi+uR*(1.0/9.0)+0.5*(uL+uR)*(t-(1.0/9.0));

    double y;
/*    if(x < x0Bar){
        y = uL;
    }else if(x > x1Bar){
        y = uR;
    }*/


    if(x < x0Bar){
        y = uL;
    }else if(x > x1Bar){
        y = uR;
    }else{
        y = uL+(uR-uL)*(x0Bar-x)/(x0Bar-x1Bar);
    }
    return y;
}


