#include "timesolver.h"
#include "thetamethod.h"

TimeSolver::TimeSolver(){

}

TimeSolver::~TimeSolver(){

}

TimeSolver* TimeSolver::Create(std::string inputFile){
    auto file = cpptoml::parse_file(inputFile);
    auto general = file->get_table("problem");
    std::string timestepping = general->get_as<std::string>("timestepping").value_or("");
    if(timestepping.compare("explicitEuler") == 0){
        return new ThetaMethod(0.0);
    }
    else{
        std::cerr<<"Invalid timesolver type"<<std::endl;
        exit(EXIT_FAILURE);
        return NULL;
    }
}
