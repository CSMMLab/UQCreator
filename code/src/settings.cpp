#include "settings.h"

Settings::Settings( std::string inputFile ) : _inputFile( inputFile ) {
    bool validConfig = true;
    try {
        _cwd      = std::filesystem::current_path();
        _inputDir = std::experimental::filesystem::path( _inputFile.parent_path() );

        auto file = cpptoml::parse_file( _inputFile );

        // section general
        auto general           = file->get_table( "general" );
        auto problemTypeString = general->get_as<std::string>( "problem" );
        if( problemTypeString ) {
            if( problemTypeString->compare( "Burgers1D" ) == 0 ) {
                _problemType = ProblemType::P_BURGERS_1D;
            }
            else if( problemTypeString->compare( "Euler1D" ) == 0 ) {
                _problemType = ProblemType::P_EULER_1D;
            }
            else if( problemTypeString->compare( "Euler2D" ) == 0 ) {
                _problemType = ProblemType::P_EULER_2D;
            }
            else {
                std::cerr << "[Inputfile][general][(string)problem] Invalid type!\nPlease set one of the following types: Burgers, Euler, Euler2D"
                          << std::endl;
                validConfig = false;
            }
        }
        else {
            std::cerr << "[Inputfile][general][(string)problem] Not set!\nPlease set one of the following types: Burgers, Euler, Euler2D"
                      << std::endl;
            validConfig = false;
        }
        auto outputDir = general->get_as<std::string>( "outputDir" );
        if( outputDir ) {
            _outputDir = *outputDir;
        }
        else {
            std::cerr << "[Inputfile][general][(string)outputDir] Not set!" << std::endl;
            validConfig = false;
        }

        // section mesh
        auto mesh          = file->get_table( "mesh" );
        auto meshDimension = mesh->get_as<unsigned>( "dimension" );
        if( meshDimension ) {
            _meshDimension = *meshDimension;
        }
        else {
            std::cerr << "[Inputfile][mesh][(unsigned)dimension] Not set!" << std::endl;
            validConfig = false;
        }
        auto outputFile = mesh->get_as<std::string>( "outputFile" );
        if( outputFile ) {
            _outputFile = _inputDir.string() + "/" + _outputDir.string() + "/" + *outputFile;
        }
        else {
            std::cerr << "[Inputfile][mesh][(string)outputFile] Not set!" << std::endl;
            validConfig = false;
        }

        // section problem
        auto problem                = file->get_table( "problem" );
        auto timesteppingTypeString = problem->get_as<std::string>( "timestepping" );
        if( timesteppingTypeString ) {
            if( timesteppingTypeString->compare( "explicitEuler" ) == 0 ) {
                _timesteppingType = TimesteppingType::T_EXPLICITEULER;
            }
            else {
                std::cerr << "[Inputfile][problem][(string)timestepping] Invalid type!\nPlease set one of the following types: explicitEuler"
                          << std::endl;
                validConfig = false;
            }
        }
        else {
            std::cerr << "[Inputfile][problem][(string)timestepping] Not set!\nPlease set one of the following types: explicitEuler" << std::endl;
            validConfig = false;
        }
        auto CFL = problem->get_as<double>( "CFL" );
        if( CFL ) {
            _CFL = *CFL;
        }
        else {
            std::cerr << "[Inputfile][problem][(double)CFL] Not set!" << std::endl;
            validConfig = false;
        }
        auto limitertypeString = problem->get_as<std::string>( "limiter" ).value_or( "none" );
        if( limitertypeString.compare( "minmod" ) == 0 ) {
            _limiterType = LimiterType::L_MINMOD;
        }
        else if( limitertypeString.compare( "none" ) == 0 ) {
            _limiterType = LimiterType::L_NONE;
        }
        else {
            std::cerr << "[Inputfile][problem][(string)limiter] Invalid type!\nPlease set one of the following types: minmod" << std::endl;
            validConfig = false;
        }
        auto tEnd = problem->get_as<double>( "tEnd" );
        if( tEnd ) {
            _tEnd = *tEnd;
        }
        else {
            std::cerr << "[Inputfile][problem][(double)tEnd] Not set!" << std::endl;
            validConfig = false;
        }

        // section moment_system
        auto moment_system     = file->get_table( "moment_system" );
        auto closureTypeString = moment_system->get_as<std::string>( "closure" );
        if( closureTypeString ) {
            if( closureTypeString->compare( "BoundedBarrier" ) == 0 ) {
                _closureType = ClosureType::C_BOUNDEDBARRIER;
            }
            else if( closureTypeString->compare( "StochasticGalerkin" ) == 0 ) {
                _closureType = ClosureType::C_STOCHASTICGALERKIN;
            }
            else if( closureTypeString->compare( "Euler" ) == 0 ) {
                _closureType = ClosureType::C_EULER_1D;
            }
            else if( closureTypeString->compare( "Euler2D" ) == 0 ) {
                _closureType = ClosureType::C_EULER_2D;
            }
            else {
                std::cerr << "[Inputfile][moment_system][(string)closure] Invalid type!\nPlease set one of the following types: BoundedBarrier, "
                             "StochasticGalerkin, Euler, Euler2D"
                          << std::endl;
                validConfig = false;
            }
        }
        else {
            std::cerr << "[Inputfile][moment_system][(string)closure] Not set!\n Please set one of the following types: BoundedBarrier, "
                         "StochasticGalerkin, Euler, Euler2D"
                      << std::endl;
            validConfig = false;
        }
        auto nMoments = moment_system->get_as<unsigned>( "moments" );
        if( nMoments ) {
            _nMoments = *nMoments;
        }
        else {
            std::cerr << "[Inputfile][moment_system][(unsigned)moments] Not set!" << std::endl;
            validConfig = false;
        }
        auto nQuadPoints = moment_system->get_as<unsigned>( "quadPoints" );
        if( nQuadPoints ) {
            _nQuadPoints = *nQuadPoints;
        }
        else {
            std::cerr << "[Inputfile][moment_system][(unsigned)quadPoints] Not set!" << std::endl;
            validConfig = false;
        }
        _maxIterations = moment_system->get_as<unsigned>( "maxIterations" ).value_or( 1000 );
        _epsilon       = moment_system->get_as<double>( "epsilon" ).value_or( 5e-5 );

        // section plot
        auto plot        = file->get_table( "plot" );
        auto plotEnabled = plot->get_as<bool>( "enabled" );
        if( plotEnabled ) {
            _plotEnabled = *plotEnabled;
        }
        else {
            std::cerr << "[Inputfile][plot][(bool)enabled] Not set!" << std::endl;
            validConfig = false;
        }
        _plotStepInterval = plot->get_as<unsigned>( "steps" ).value_or( 0 );
        _plotTimeInterval = plot->get_as<double>( "time" ).value_or( -1.0 );

    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _inputFile << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }
    if( !validConfig ) {
        exit( EXIT_FAILURE );
    }
}

Settings::~Settings() {}

// general
ProblemType Settings::GetProblemType() const { return _problemType; }
unsigned Settings::GetNStates() const { return _nStates; }
void Settings::SetNStates( unsigned n ) { _nStates = n; }
std::string Settings::GetInputFile() const { return _inputFile; }
std::string Settings::GetInputDir() const { return _inputDir; }
std::string Settings::GetOutputDir() const { return _outputDir; }

// mesh
unsigned Settings::GetMeshDimension() const { return _meshDimension; }
unsigned Settings::GetNumCells() const { return _numCells; }
void Settings::SetNumCells( unsigned n ) { _numCells = n; }
std::string Settings::GetOutputFile() const { return _outputFile; }

// problem
TimesteppingType Settings::GetTimesteppingType() const { return _timesteppingType; }
double Settings::GetCFL() const { return _CFL; }
double Settings::GetTEnd() const { return _tEnd; }
double Settings::GetGamma() const { return _gamma; }
void Settings::SetGamma( double gamma ) { _gamma = gamma; }

// moment_system
ClosureType Settings::GetClosureType() const { return _closureType; }
unsigned Settings::GetNMoments() const { return _nMoments; }
unsigned Settings::GetNQuadPoints() const { return _nQuadPoints; }
LimiterType Settings::GetLimiterType() const { return _limiterType; }
unsigned Settings::GetMaxIterations() const { return _maxIterations; }
double Settings::GetEpsilon() const { return _epsilon; }

// plot
bool Settings::GetPlotEnabled() const { return _plotEnabled; }
unsigned Settings::GetPlotStepInterval() const { return _plotStepInterval; }
double Settings::GetPlotTimeInterval() const { return _plotTimeInterval; }
