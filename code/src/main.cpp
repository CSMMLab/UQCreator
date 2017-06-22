#include <blaze/Blaze.h>
#include <cpptoml.h>
#include "burgerssolver.h"
#include <omp.h>

int main(int argc, char* argv[]){
    std::string configFile = "";
/*
    std::string usage_help =
        "\n"
        "Usage: " + std::string(argv[0]) + " -i inputfile\n\n"
        "Options:\n"
        "  -h               displays this message\n"
        "  -i <config.json> provide toml input file <input.toml>\n"
        "  -t               number of threads to use (default = all available)\n";

    if (argc < 3)
    {
        std::cout << usage_help;
        return 0;
    }
    for (int i=1; i<argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-h")
        {
            std::cout << usage_help;
            //return 0;
        }
        else if (arg == "-i")
        {
            configFile = std::string(argv[++i]);
        }
        else if (arg == "-t")
        {
             omp_set_num_threads(std::atoi(argv[++i]));
        }
        else{
            //std::cout << usage_help;
            //return 0;
        }
    }*/

    int nCells = 100;
    double tEnd = 0.1;
    double cfl = 0.8;
    double a = 0.0;
    double b = 3.0;
    double uL = 12.0;
    double uR = 3.0;
    BurgersSolver s(nCells, tEnd, cfl, a, b,uL,uR);
    s.Solve();
    s.Print();


    return 0;
}
