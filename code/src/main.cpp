#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <iostream>

#include "problem.h"
#include "burgers.h"
#include "momentsolver.h"

bool CheckInput(std::string& configFile, int argc, char* argv[]){
    std::string usage_help =
        "\n"
        "Usage: " + std::string(argv[0]) + " -c inputfile\n\n"
        "Options:\n"
        "  -h               displays this message\n"
        "  -c <config.json> provide toml input file <input.toml>\n"
        "  -t               number of threads to use (default = all available)\n";

    if (argc < 3){
        std::cout << usage_help;
        return false;
    }
    for (int i=1; i<argc; i++){
        std::string arg = argv[i];
        if (arg == "-h"){
            std::cout << usage_help;
            return false;
        }
        else if (arg == "-c"){
            configFile = std::string(argv[++i]);
            std::ifstream f(configFile);
            if(!f.is_open())
                return false;
        }
        else if (arg == "-t"){
             omp_set_num_threads(std::atoi(argv[++i]));
        }
        else{
            std::cout << usage_help;
            return false;
        }
    }
    return true;
}

void PrintInit(std::string configFile){
    std::cout << "UQCreator" << std::endl;
    std::cout << "==================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Config file:\t" + configFile << std::endl;
    std::cout << "==================================" << std::endl;
    std::ifstream ifs(configFile);
    if(ifs.is_open()){
        std::string line;
        while(!ifs.eof()){
            std::getline(ifs,line);
            std::cout << " " << line << std::endl;
        }
    }
    std::cout << "==================================\n" << std::endl;
}

int main(int argc, char* argv[]){
    std::string configFile = "";

    if(!CheckInput(configFile, argc, argv)){return -1;}
    PrintInit(configFile);

    Problem* problem = Problem::Create(configFile);
    MomentSolver* solver = new MomentSolver(problem);

    solver->Solve();
    solver->Plot();

    delete solver;
    delete problem;

    std::cout << "\nProcess exited normally." << std::endl;

    return 0;
}
