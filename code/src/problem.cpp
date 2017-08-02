#include "problem.h"
#include "burgers.h"
#include "euler.h"

Problem::Problem(std::string inputFile) : _inputFile(inputFile){
    try{
        auto file = cpptoml::parse_file(_inputFile);

        auto general = file->get_table("general");
        _outputFolder = general->get_as<std::string>("outputFolder").value_or("");

        _mesh = new Mesh(_inputFile);

        auto problem = file->get_table("problem");
        _CFL = problem->get_as<double>("CFL").value_or(-1.0);
        _limiter = problem->get_as<std::string>("limiter").value_or("none");
        _tEnd = problem->get_as<double>("tEnd").value_or(-1.0);

        auto momentSystem = file->get_table("moment system");
        std::string quadType = momentSystem->get_as<std::string>("quadType").value_or("none");
        if(quadType.compare("legendre"))
            _quadType = QUAD_TYPE_LEGENDRE;
        else if(quadType.compare("hermite"))
            _quadType = QUAD_TYPE_HERMITE;
        _nQuadPoints = momentSystem->get_as<int>("quadPoints").value_or(-1);
        _nMoments = momentSystem->get_as<int>("moments").value_or(-1);
        _maxIterations = momentSystem->get_as<int>("maxIterations").value_or(-1);
        _epsilon = momentSystem->get_as<double>("epsilon").value_or(-1.0);\
    }
    catch (const cpptoml::parse_exception& e){
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

Problem::~Problem(){

}

Problem* Problem::Create(std::string inputFile){
    auto file = cpptoml::parse_file(inputFile);
    auto general = file->get_table("general");
    std::string problem = general->get_as<std::string>("problem").value_or("");
    if(problem.compare("Burgers") == 0){
        return new Burgers(inputFile);
    }
    else{
        std::cerr<<"Invalid problem type"<<std::endl;
        exit(EXIT_FAILURE);
        return NULL;
    }
}

int Problem::GetQuadType(){
    return _quadType;
}

int Problem::GetNQuadPoints(){
    return _nQuadPoints;
}

int Problem::GetNMoments(){
    return _nMoments;
}

int Problem::GetMaxIterations(){
    return _maxIterations;
}

double Problem::GetEpsilon(){
    return _epsilon;
}

double Problem::GetCFL(){
    return _CFL;
}

double Problem::GetTEnd(){
    return _tEnd;
}

Mesh* Problem::GetMesh(){
    return _mesh;
}
