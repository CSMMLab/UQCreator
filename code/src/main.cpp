#include <QApplication>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <iostream>

#include "mesh.h"
#include "momentsolver.h"
#include "problem.h"
#include "settings.h"

bool CheckInput( std::string& configFile, bool& batchMode, int argc, char* argv[] ) {
    std::string usage_help = "\n"
                             "Usage: " +
                             std::string( argv[0] ) +
                             " -c inputfile\n\n"
                             "Options:\n"
                             "  -h               displays this message\n"
                             "  -c <config.json> provide toml input file <input.toml>\n"
                             "  -b               runs in batch mode (not displaying plots)";

    if( argc < 3 ) {
        std::cout << usage_help;
        return false;
    }
    for( int i = 1; i < argc; i++ ) {
        std::string arg = argv[i];
        if( arg == "-h" ) {
            std::cout << usage_help;
            return false;
        }
        else if( arg == "-c" ) {
            configFile = std::string( argv[++i] );
            std::ifstream f( configFile );
            if( !f.is_open() ) return false;
        }
        else if( arg == "-b" ) {
            batchMode = true;
        }
        else {
            std::cout << usage_help;
            return false;
        }
    }
    return true;
}

void PrintInit( std::string configFile ) {
    std::cout << "UQCreator" << std::endl;
    std::cout << "==================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Config file:\t" + configFile << std::endl;
    std::cout << "==================================" << std::endl;
    std::ifstream ifs( configFile );
    if( ifs.is_open() ) {
        std::string line;
        while( !ifs.eof() ) {
            std::getline( ifs, line );
            std::cout << " " << line << std::endl;
        }
    }
    std::cout << "==================================\n" << std::endl;
}

int main( int argc, char* argv[] ) {
    std::string configFile = "";
    bool batchMode         = false;

    if( !CheckInput( configFile, batchMode, argc, argv ) ) {
        return -1;
    }
    PrintInit( configFile );

    QApplication* app = nullptr;
    if( !batchMode ) {
        app = new QApplication( argc, argv );
    }

    Settings* settings   = new Settings( configFile );
    Mesh* mesh           = Mesh::Create( settings );
    Problem* problem     = Problem::Create( settings );
    MomentSolver* solver = new MomentSolver( settings, mesh, problem );

    solver->Solve();

    std::cout << "\nProcess exited normally." << std::endl;

    delete solver;
    delete problem;
    delete mesh;
    delete settings;

    if( !batchMode && app ) {
        return app->exec();
    }
    else {
        return EXIT_SUCCESS;
    }
}
