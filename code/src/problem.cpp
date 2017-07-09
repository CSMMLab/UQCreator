#include "problem.h"
#include "burgers.h"
#include "euler.h"

Problem::Problem(std::string inputFile) : _inputFile(inputFile){
    try{
        auto file = cpptoml::parse_file(_inputFile);

        auto general = file->get_table("general");
        _outputFolder = general->get_as<std::string>("outputFolder").value_or("");

        _mesh = new Mesh(_inputFile);
    }
    catch (const cpptoml::parse_exception& e){
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

Problem::~Problem(){

}

Problem* Problem::create(std::string inputFile){
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

int Problem::getQuadType(){
    return _quadType;
}

int Problem::getNQuadPoints(){
    return _nQuadPoints;
}
