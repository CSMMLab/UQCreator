#include "burgers.h"

Burgers::Burgers(std::string inputFile) : Problem(inputFile)
{
    try{
        auto file = cpptoml::parse_file(_inputFile);

        auto problem = file->get_table("burgers");
        _CFL = problem->get_as<double>("CFL").value_or(-1.0);
        _limiter = problem->get_as<std::string>("limiter").value_or("none");
        _tEnd = problem->get_as<double>("tEnd").value_or(-1.0);
    }
    catch (const cpptoml::parse_exception& e){
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    /*
    _u = new double[_nCells+4];
    _x = _mesh->getGrid();
    _dx = _mesh->getSpacing();
    _dt = _timeDiscretization->getDt(0.8);
    _nCells = _mesh->getNumCells();
    _nTimeSteps = tEnd/_dt;
    for( int j = 0; j<_nCells+4; ++j){
        _u[j] = IC(_x[j],uL,uR);
    }
    */
}

double Burgers::H(double u, double v, double w){
    return v-(_dt/_dx)*(g(v,w)-g(u,v));
}

double Burgers::g(double u, double v){
    return f(u);
}

double Burgers::f(double u){
    return 0.5*u*u;
}

void Burgers::solve(){
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

void Burgers::print() const{
    std::ofstream out("outFile");
    for( int j = 2; j < _nCells+2; ++j){
        out<<_x[j]<<" "<<_u[j]<<std::endl;
    }
}

void Burgers::plot() const{
    /*
    std::vector<double> x,u;
    x.assign(_x, _x + _nCells+4);
    u.assign(_u, _u + _nCells+4);
    Gnuplot gp;
    gp << "plot '-' with lines notitle\n";
    gp << "set xlabel 'Space'\n";
    gp << "set ylabel 'Velocity\n";
    gp.send1d(std::make_pair(x, u));
    */
}

void Burgers::writeToFile(std::string filename, int filetype) const{

}


