#include "settings.h"

#include <mpi.h>

Settings::Settings( std::string inputFile ) : _inputFile( inputFile ), _numDimXi( 1 ) {
    auto log = spdlog::get( "event" );

    bool validConfig = true;
    try {
        int ierr;
        ierr      = MPI_Comm_rank( MPI_COMM_WORLD, &_mype );
        ierr      = MPI_Comm_size( MPI_COMM_WORLD, &_npes );
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
                log->error( "[inputfile] [general] 'problem' is invalid!\nPlease set one of the following types: Burgers, Euler, Euler2D" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [general] 'problem' not set!\nPlease set one of the following types: Burgers, Euler, Euler2D" );
            validConfig = false;
        }
        auto outputDir = general->get_as<std::string>( "outputDir" );
        if( outputDir ) {
            _outputDir = *outputDir;
        }
        else {
            log->error( "[inputfile] [general] 'outputDir' not set!" );
            validConfig = false;
        }

        // section mesh
        auto mesh          = file->get_table( "mesh" );
        auto meshDimension = mesh->get_as<unsigned>( "dimension" );
        if( meshDimension ) {
            _meshDimension = *meshDimension;
        }
        else {
            log->error( "[inputfile] [mesh] 'dimension' not set!" );
            validConfig = false;
        }
        auto outputFile = mesh->get_as<std::string>( "outputFile" );
        if( outputFile ) {
            _outputFile = _inputDir.string() + "/" + _outputDir.string() + "/" + *outputFile;
        }
        else {
            log->error( "[inputfile] [mesh] 'outputFile' not set!" );
            validConfig = false;
        }

        auto continueFile = mesh->get_as<std::string>( "continueFile" );
        if( continueFile ) {
            _continueFile = _inputDir.string() + "/" + *continueFile;
        }

        // section problem
        auto problem                = file->get_table( "problem" );
        auto timesteppingTypeString = problem->get_as<std::string>( "timestepping" );
        if( timesteppingTypeString ) {
            if( timesteppingTypeString->compare( "explicitEuler" ) == 0 ) {
                _timesteppingType = TimesteppingType::T_EXPLICITEULER;
            }
            else {
                log->error( "[inputfile] [problem] 'timestepping' is invalid!\nPlease set one of the following types: explicitEuler" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [problem] 'timestepping' not set!\nPlease set one of the following types: explicitEuler" );
            validConfig = false;
        }
        auto distribution = problem->get_as<std::string>( "distribution" );
        if( distribution ) {
            if( distribution->compare( "Legendre" ) == 0 ) {
                _distributionType = DistributionType::D_LEGENDRE;
            }
            else if( distribution->compare( "Hermite" ) == 0 ) {
                _distributionType = DistributionType::D_HERMITE;
            }
            else {
                log->error( "[inputfile] [problem] 'distribution' is invalid!\nPlease set one of the following types: Legendre, Hermite" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [problem] 'timestepping' not set!\nPlease set one of the following types: explicitEuler" );
            validConfig = false;
        }
        auto CFL = problem->get_as<double>( "CFL" );
        if( CFL ) {
            _CFL = *CFL;
        }
        else {
            log->error( "[inputfile] [problem] 'CFL' not set!" );
            validConfig = false;
        }
        auto tEnd = problem->get_as<double>( "tEnd" );
        if( tEnd ) {
            _tEnd = *tEnd;
        }
        else {
            log->error( "[inputfile] [problem] 'tEnd' not set!" );
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
                log->error( "[inputfile] [moment_system] 'closure' is invalid!\nPlease set one of the following types: BoundedBarrier, "
                            "StochasticGalerkin, Euler, Euler2D" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [moment_system] 'closure' not set!\n Please set one of the following types: BoundedBarrier, "
                        "StochasticGalerkin, Euler, Euler2D" );
            validConfig = false;
        }
        auto nMoments = moment_system->get_as<unsigned>( "moments" );
        if( nMoments ) {
            _nMoments = *nMoments;
        }
        else {
            log->error( "[inputfile] [moment_system] 'moments' not set!" );
            validConfig = false;
        }
        auto nQuadPoints = moment_system->get_as<unsigned>( "quadPoints" );
        if( nQuadPoints ) {
            _nQuadPoints = *nQuadPoints;
            _nQTotal     = unsigned( std::pow( _nQuadPoints, _numDimXi ) );
            // determine size of quad array on PE
            int nQPE = int( ( int( _nQTotal ) - 1 ) / _npes ) + 1;
            if( _mype == _npes - 1 ) {
                nQPE = int( _nQTotal ) - _mype * nQPE;
                if( nQPE < 0 ) {
                    nQPE = 0;
                }
            }
            _nQPE   = unsigned( nQPE );
            _kStart = _mype * ( ( _nQTotal - 1 ) / _npes + 1.0 );
            _kEnd   = _kStart + _nQPE - 1;
        }
        else {
            log->error( "[inputfile] [moment_system] 'quadPoints' not set!" );
            validConfig = false;
        }
        _maxIterations = moment_system->get_as<unsigned>( "maxIterations" ).value_or( 1000 );
        _epsilon       = moment_system->get_as<double>( "epsilon" ).value_or( 5e-5 );

    } catch( const cpptoml::parse_exception& e ) {
        log->error( "Failed to parse {0}: {1}", _inputFile.c_str(), e.what() );
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
bool Settings::HasContinueFile() const { return !_continueFile.empty(); }
std::string Settings::GetContinueFile() const { return _continueFile; }

// problem
TimesteppingType Settings::GetTimesteppingType() const { return _timesteppingType; }
double Settings::GetCFL() const { return _CFL; }
double Settings::GetTEnd() const { return _tEnd; }
unsigned Settings::GetNDimXi() const { return _numDimXi; }
double Settings::GetGamma() const { return _gamma; }
void Settings::SetGamma( double gamma ) { _gamma = gamma; }
DistributionType Settings::GetDistributionType() const { return _distributionType; }

// moment_system
ClosureType Settings::GetClosureType() const { return _closureType; }
unsigned Settings::GetNMoments() const { return _nMoments; }
unsigned Settings::GetNQuadPoints() const { return _nQuadPoints; }
unsigned Settings::GetNQTotal() const { return _nQTotal; }
LimiterType Settings::GetLimiterType() const { return _limiterType; }
unsigned Settings::GetMaxIterations() const { return _maxIterations; }
void Settings::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }
double Settings::GetEpsilon() const { return _epsilon; }

// plot
unsigned Settings::GetPlotStepInterval() const { return _plotStepInterval; }
double Settings::GetPlotTimeInterval() const { return _plotTimeInterval; }

// MPI
int Settings::GetMyPE() const { return _mype; }
int Settings::GetNPEs() const { return _npes; }
unsigned Settings::GetKStart() const { return _kStart; }
unsigned Settings::GetKEnd() const { return _kEnd; }
unsigned Settings::GetNqPE() const { return _nQPE; }
