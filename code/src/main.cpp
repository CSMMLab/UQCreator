#define LOG(x) do { std::clog << x << std::endl; } while (0)
#define DEBUG(x) do { std::cerr << x << std::endl; } while (0)

#include "problem.h"

int main(int argc, char* argv[]){
    std::string configFile = "";

    std::string usage_help =
        "\n"
        "Usage: " + std::string(argv[0]) + " -c inputfile\n\n"
        "Options:\n"
        "  -h               displays this message\n"
        "  -c <config.json> provide toml input file <input.toml>\n"
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
            return 0;
        }
        else if (arg == "-c")
        {
            configFile = std::string(argv[++i]);
        }
        else if (arg == "-t")
        {
             omp_set_num_threads(std::atoi(argv[++i]));
        }
        else{
            std::cout << usage_help;
            return 0;
        }
    }

    LOG("UQCreator");
    LOG("================================\n");
    LOG("Using config file \t: " + configFile);

    Problem* problem = Problem::Create(configFile);

    delete problem;

    LOG("\nProcess exits normally.");


    return 0;
}
