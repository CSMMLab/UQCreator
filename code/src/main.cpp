#include <blaze/Blaze.h>
#include <cpptoml.h>

int main(int argc, char* argv[]){
    std::string inputFile = "";

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
            return 0;
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
            std::cout << usage_help;
            return 0;
        }
    }

    return 0;
}
