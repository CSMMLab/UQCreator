#include "problem.h"
#include "burgers.h"

Problem::Problem( std::string inputFile ) : _inputFile( inputFile ) {
    try {
        auto file = cpptoml::parse_file( _inputFile );

        auto general = file->get_table( "general" );
        _outputDir   = general->get_as<std::string>( "outputDir" ).value_or( "" );

        _mesh = new Mesh( _inputFile );

        auto problem = file->get_table( "problem" );
        _CFL         = problem->get_as<double>( "CFL" ).value_or( -1.0 );
        _limiter     = problem->get_as<std::string>( "limiter" ).value_or( "none" );
        _tEnd        = problem->get_as<double>( "tEnd" ).value_or( -1.0 );
        _nStates     = problem->get_as<double>( "nStates" ).value_or( -1.0 );

        auto momentSystem    = file->get_table( "moment_system" );
        std::string quadType = momentSystem->get_as<std::string>( "quadType" ).value_or( "none" );
        if( quadType.compare( "legendre" ) )
            _quadType = QUAD_TYPE_LEGENDRE;
        else if( quadType.compare( "hermite" ) )
            _quadType = QUAD_TYPE_HERMITE;
        _nQuadPoints   = momentSystem->get_as<int>( "quadPoints" ).value_or( -1 );
        _nMoments      = momentSystem->get_as<int>( "moments" ).value_or( -1 );
        _maxIterations = momentSystem->get_as<int>( "maxIterations" ).value_or( -1 );
        _epsilon       = momentSystem->get_as<double>( "epsilon" ).value_or( -1.0 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
}

Problem::~Problem() { delete _mesh; }

Problem* Problem::Create( std::string inputFile ) {
    auto file           = cpptoml::parse_file( inputFile );
    auto general        = file->get_table( "general" );
    std::string problem = general->get_as<std::string>( "problem" ).value_or( "" );
    if( problem.compare( "Burgers" ) == 0 ) {
        return new Burgers( inputFile );
    }
    else {
        std::cerr << "Invalid problem type" << std::endl;
        exit( EXIT_FAILURE );
        return NULL;
    }
}

int Problem::GetQuadType() const { return _quadType; }

int Problem::GetNQuadPoints() const { return _nQuadPoints; }

int Problem::GetNMoments() const { return _nMoments; }

int Problem::GetNStates() const { return _nStates; }

int Problem::GetMaxIterations() const { return _maxIterations; }

double Problem::GetEpsilon() const { return _epsilon; }

double Problem::GetCFL() const { return _CFL; }

double Problem::GetTEnd() const { return _tEnd; }

Mesh* Problem::GetMesh() const { return _mesh; }

std::string Problem::GetInputFile() const { return _inputFile; }

std::string Problem::GetLimiter() const { return _limiter; }

std::string Problem::GetOutputDir() const { return _outputDir; }
