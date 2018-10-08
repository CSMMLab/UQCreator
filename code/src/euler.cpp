#include "euler.h"

Euler::Euler( std::string inputFile ) : Problem( inputFile ) {
    try {
        auto file = cpptoml::parse_file( _inputFile );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

double Euler::GetGamma() const { return _gamma; }

void Euler::Solve() {}

void Euler::Print() const {}

void Euler::Plot() const {}

void Euler::WriteToFile( std::string filename, int filetype ) const {}
