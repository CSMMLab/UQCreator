#include "problem.h"
#include "burgers.h"
#include "euler.h"
#include "euler2d.h"

Problem::Problem( const Settings* settings ) : _settings( settings ) {}

Problem::~Problem() {}

Problem* Problem::Create( const Settings* settings ) {
    auto file           = cpptoml::parse_file( settings->GetInputFile() );
    auto general        = file->get_table( "general" );
    std::string problem = general->get_as<std::string>( "problem" ).value_or( "" );
    if( problem.compare( "Burgers" ) == 0 ) {
        return new Burgers( settings );
    }
    else if( problem.compare( "Euler" ) == 0 ) {
        return new Euler( settings );
    }
    else if( problem.compare( "Euler2D" ) == 0 ) {
        return new Euler2D( settings );
    }
    else {
        std::cerr << "Invalid problem type!" << std::endl;
        exit( EXIT_FAILURE );
    }
}
