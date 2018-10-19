#include "settings.h"

Settings::Settings( std::string inputFile ) {
    try {
        _inputDir = std::experimental::filesystem::path( inputFile );
        auto file = cpptoml::parse_file( _inputFile );

        auto general           = file->get_table( "general" );
        auto problemTypeString = general->get_as<std::string>( "problem" ).value_or( "none" );
        if( problemTypeString.compare( "Burgers1D" ) == 0 ) {
            _problemType = ProblemType::P_BURGERS_1D;
        }
        else if( problemTypeString.compare( "Euler1D" ) == 0 ) {
            _problemType = ProblemType::P_EULER_1D;
        }
        else if( problemTypeString.compare( "Euler2D" ) == 0 ) {
            _problemType = ProblemType::P_EULER_2D;
        }
        else {
            std::cerr << "Invalid problem type!" << std::endl;
        }

        auto problem           = file->get_table( "problem" );
        _CFL                   = problem->get_as<double>( "CFL" ).value_or( -1.0 );
        auto limitertypeString = problem->get_as<std::string>( "limiter" ).value_or( "none" );
        if( limitertypeString.compare( "minmod" ) == 0 ) {
            _limiter = LimiterType::L_MINMOD;
        }
        else if( limitertypeString.compare( "none" ) == 0 ) {
            _limiter = LimiterType::L_NONE;
        }
        else {
            std::cerr << "Invalid limiter type!" << std::endl;
        }

        _tEnd = problem->get_as<double>( "tEnd" ).value_or( -1.0 );

    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

Settings::~Settings() {}

unsigned Settings::GetNQuadPoints() const { return _nQuadPoints; }

unsigned Settings::GetNMoments() const { return _nMoments; }

unsigned Settings::GetNStates() const { return _nStates; }

unsigned Settings::GetMaxIterations() const { return _maxIterations; }

double Settings::GetEpsilon() const { return _epsilon; }

double Settings::GetCFL() const { return _CFL; }

double Settings::GetTEnd() const { return _tEnd; }

std::string Settings::GetInputFile() const { return _inputFile; }

LimiterType Settings::GetLimiter() const { return _limiter; }

std::string Settings::GetOutputDir() const { return _outputDir; }

ClosureType Settings::GetClosureType() const { return _closureType; }

ProblemType Settings::GetProblemType() const { return _problemType; }

unsigned Settings::GetNumCells() const { return _numCells; }

void Settings::SetNumCells( unsigned n ) { _numCells = n; }

double Settings::GetGamma() const { return _gamma; }
