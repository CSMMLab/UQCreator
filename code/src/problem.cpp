#include "problem.h"

Problem::Problem(std::string inputFile) : _inputFile(inputFile){

}

Problem::~Problem(){
    delete _discretization;
}

Problem* Problem::create(std::string inputFile){
    auto file = cpptoml::parse_file(inputFile);
    auto general = file->get_table("general");
    std::string problem = general->get_as<std::string>("problem").value_or("");
    /*
    if(problem.compare("Burgers") == 0)
        return new Burgers(inputFile);
    else{
        std::cerr<<"Invalid problem type"<<std::endl;
        exit(EXIT_FAILURE);
        return NULL;
    }
    */
}

void Problem::parse(){
    try{
        auto file = cpptoml::parse_file(_inputFile);

        auto general = file->get_table("general");
        _outputFolder = general->get_as<std::string>("outputFolder").value_or("");

        auto spaceTime = file->get_table("space-time");
        _dimension = spaceTime->get_as<int>("dimension").value_or(-1);
        _discretization = new int[_dimension];
        switch(_dimension){
            case 1: _discretization[0] = spaceTime->get_as<int>("xCells").value_or(-1);
            case 2: _discretization[1] = spaceTime->get_as<int>("yCells").value_or(-1);
            case 3: _discretization[2] = spaceTime->get_as<int>("zCells").value_or(-1);
            default: cpptoml::parse_exception("Invalid dimension");
        }
        _tEnd = spaceTime->get_as<double>("tEnd").value_or(-1.0);
        _CFL = spaceTime->get_as<double>("CFL").value_or(-1.0);
        _limiter = spaceTime->get_as<std::string>("limiter").value_or("none");

    }
    catch (const cpptoml::parse_exception& e){
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}


