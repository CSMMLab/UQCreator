#include <mpi.h>

#include "settings.h"

namespace cpptoml {
inline std::shared_ptr<table> parse_file( const std::istringstream& content ) {
    parser p{const_cast<std::istringstream&>( content )};
    return p.parse();
}
}    // namespace cpptoml

Settings::Settings( std::string inputFile ) : _inputFile( inputFile ), _hasExactSolution( false ) {
    auto file = cpptoml::parse_file( _inputFile );
    Init( file, false );
}

Settings::Settings( const std::istringstream& inputStream ) {
    auto file = cpptoml::parse_file( inputStream );
    Init( file, true );
}

void Settings::Init( std::shared_ptr<cpptoml::table> file, bool restart ) {
    auto log = spdlog::get( "event" );

    bool validConfig = true;
    try {
        MPI_Comm_rank( MPI_COMM_WORLD, &_mype );
        MPI_Comm_size( MPI_COMM_WORLD, &_npes );
        _cwd      = std::filesystem::current_path();
        _inputDir = std::experimental::filesystem::path( _inputFile.parent_path() );

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
            else if( problemTypeString->compare( "ShallowWater1D" ) == 0 ) {
                _problemType = ProblemType::P_SHALLOWWATER_1D;
            }
            else if( problemTypeString->compare( "ShallowWater2D" ) == 0 ) {
                _problemType = ProblemType::P_SHALLOWWATER_2D;
            }
            else if( problemTypeString->compare( "PNEquations2D" ) == 0 ) {
                _problemType = ProblemType::P_PNEQUATIONS_2D;
            }
            else {
                log->error( "[inputfile] [general] 'problem' is invalid!\nPlease set one of the following types: Burgers1D, Euler1D, "
                            "Euler2D,ShallowWater1D,ShallowWater2D" );
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

        auto restartFile = general->get_as<std::string>( "restartFile" );
        if( restartFile ) {
            _restartFile = _inputDir.string() + "/" + *restartFile;
            if( restart ) _loadLambda = general->get_as<bool>( "importDualState" ).value_or( false );
        }
        if( !restart ) _loadLambda = general->get_as<bool>( "importDualState" ).value_or( false );

        auto icFile = general->get_as<std::string>( "icFile" );
        if( icFile ) {
            _icFile = _inputDir.string() + "/" + *icFile;
        }

        if( !restart ) {
            auto refFile    = general->get_as<std::string>( "referenceSolution" );
            _writeFrequency = general->get_as<int64_t>( "writeFrequency" ).value_or( 1000 );
            if( refFile ) {
                _referenceFile = _inputDir.string() + "/" + *refFile;
            }
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
        auto distribution = problem->get_array_of<cpptoml::array>( "distribution" );
        if( distribution ) {
            _numDimXi = unsigned( distribution->at( 0 )->get_array_of<std::string>()->size() );
            _distributionType.resize( _numDimXi );
            _sigma.resize( _numDimXi );
            auto dist  = ( *distribution )[0]->get_array_of<std::string>();
            auto sigma = ( *distribution )[1]->get_array_of<double>();
            for( unsigned i = 0; i < dist->size(); ++i ) {
                if( dist->at( i ).compare( "Legendre" ) == 0 ) {
                    _distributionType[i] = DistributionType::D_LEGENDRE;
                }
                else if( dist->at( i ).compare( "Hermite" ) == 0 ) {
                    _distributionType[i] = DistributionType::D_HERMITE;
                }
                else {
                    log->error( "[inputfile] [problem] 'distribution' is invalid!\nPlease set one of the following types: Legendre, Hermite" );
                    validConfig = false;
                }
                _sigma[i] = ( *sigma )[i];
            }
        }
        else {
            log->error( "[inputfile] [problem] 'distribution' not set!\nPlease set one of the following types: Legendre, Hermite" );
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
        _minResidual = problem->get_as<double>( "residual" ).value_or( -1.0 );

        // section moment_system
        auto moment_system     = file->get_table( "moment_system" );
        auto closureTypeString = moment_system->get_as<std::string>( "closure" );
        if( closureTypeString ) {
            if( closureTypeString->compare( "BoundedBarrier" ) == 0 ) {
                _closureType = ClosureType::C_BOUNDEDBARRIER;
            }
            else if( closureTypeString->compare( "LogSin" ) == 0 ) {
                _closureType = ClosureType::C_LOGSIN;
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
            else if( closureTypeString->compare( "ShallowWater" ) == 0 ) {
                _closureType = ClosureType::C_SHALLOWWATER_1D;
            }
            else if( closureTypeString->compare( "ShallowWater2D" ) == 0 ) {
                _closureType = ClosureType::C_SHALLOWWATER_2D;
            }
            else if( closureTypeString->compare( "L2Filter" ) == 0 ) {
                _closureType = ClosureType::C_L2FILTER;
            }
            else if( closureTypeString->compare( "LassoFilter" ) == 0 ) {
                _closureType = ClosureType::C_LASSOFILTER;
            }
            else {
                log->error( "[inputfile] [moment_system] 'closure' is invalid!\nPlease set one of the following types: BoundedBarrier, LogSin, "
                            "StochasticGalerkin, Euler, Euler2D,L2Filter,LassoFilter" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [moment_system] 'closure' not set!\n Please set one of the following types: BoundedBarrier, LogSin, "
                        "StochasticGalerkin, Euler, Euler2D,L2Filter,LassoFilter" );
            validConfig = false;
        }
        auto momentSettings = moment_system->get_array_of<cpptoml::array>( "moments" );

        if( momentSettings ) {
            auto momentArray       = ( *momentSettings )[1]->get_array_of<int64_t>();
            auto degreeType        = ( *momentSettings )[0]->get_array_of<std::string>();
            _nRefinementLevels     = unsigned( momentArray->size() );
            _nTotalRefinementLevel = VectorU( _nRefinementLevels );
            _refinementLevel       = VectorU( _nRefinementLevels );
            for( unsigned i = 0; i < _nRefinementLevels; ++i ) _refinementLevel[i] = unsigned( ( *momentArray )[i] );
            _nMoments = unsigned( ( *momentArray )[momentArray->size() - 1] );    // unsigned( ( *nMoments )[nMoments->size() - 1] );
            // compute nTotal
            if( degreeType->at( 0 ).compare( "maxDegree" ) == 0 ) {
                _useMaxDegree = true;
            }
            else if( degreeType->at( 0 ).compare( "totalDegree" ) == 0 ) {
                _useMaxDegree = false;
            }
            else {
                log->error( "[inputfile] [moment_system] 'moments' is invalid!\nPlease set one of the following types: totalDegree, maxDegree" );
                validConfig = false;
            }
            _nTotal = 0;

            // setup map from k (0,...,_nQTotal-1) to individual indices
            std::vector<std::vector<unsigned>> _polyIndices;

            // setup map from i\in(0,...,nTotal-1) to individual indices for basis function calculation
            VectorU nTotal( _nRefinementLevels );
            std::vector<unsigned> indexTest;
            indexTest.resize( _numDimXi );
            int totalDegree;    // compute total degree of basis function i
            int previousDegree = -1;
            // loop over all levels and only store indices of certain level to ensure correct ordering
            for( unsigned level = 0; level < _nRefinementLevels; ++level ) {
                for( unsigned i = 0; i < std::pow( _nMoments + 1, _numDimXi ); ++i ) {
                    totalDegree = 0;
                    for( unsigned l = 0; l < _numDimXi; ++l ) {
                        indexTest[l] = unsigned( ( i - i % unsigned( std::pow( _nMoments + 1, l ) ) ) / unsigned( std::pow( _nMoments + 1, l ) ) ) %
                                       ( _nMoments + 1 );
                        totalDegree += indexTest[l];
                    }
                    // if total degree is sufficiently small or max degree is used, indices are stored
                    if( ( unsigned( totalDegree ) <= GetPolyDegreeforRefLevel( level ) && totalDegree > previousDegree ) || this->UsesMaxDegree() ) {
                        _polyIndices.push_back( indexTest );
                        unsigned ii = _polyIndices.size() - 1;
                        for( unsigned j = 0; j < _polyIndices[ii].size(); ++j ) {
                            std::cout << _polyIndices[ii][j] << " ";
                        }
                        std::cout << ", degree " << totalDegree << std::endl;
                    }
                }
                if( UsesMaxDegree() ) break;
                previousDegree                = int( GetPolyDegreeforRefLevel( level ) );
                _nTotalRefinementLevel[level] = _polyIndices.size();
            }
        }
        else {
            log->error( "[inputfile] [moment_system] 'moments' not set!" );
            validConfig = false;
        }

        auto quadratureSettings = moment_system->get_array_of<cpptoml::array>( "quadPoints" );

        if( quadratureSettings ) {
            auto nQArray        = ( *quadratureSettings )[1]->get_array_of<int64_t>();
            auto quadratureType = ( *quadratureSettings )[0]->get_array_of<std::string>();
            _nQuadPoints        = unsigned( ( *nQArray )[nQArray->size() - 1] );

            // computing Nq Total
            if( quadratureType->at( 0 ).compare( "sparseGrid" ) == 0 ) {
                _gridType = G_SPARSEGRID;
                _nQTotal  = unsigned( std::pow( 2, _nQuadPoints ) ) + 1u;
                _nQTotal  = 13;
            }
            else {    // tensorizedGrid
                _gridType = G_TENSORIZEDGRID;
                _nQTotal  = unsigned( std::pow( _nQuadPoints, _numDimXi ) );
            }
            /*
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
            _kEnd   = _kStart + _nQPE - 1;*/
            _kEnd   = unsigned( std::floor( ( double( _mype ) + 1.0 ) * double( _nQTotal / double( _npes ) ) ) );
            _kStart = unsigned( std::ceil( double( _mype ) * ( double( _nQTotal ) / double( _npes ) ) ) );
            if( unsigned( std::ceil( ( double( _mype ) + 1.0 ) * ( double( _nQTotal ) / double( _npes ) ) ) ) == _kEnd ) _kEnd = _kEnd - 1;
            _nQPE = _kEnd - _kStart + 1;
        }
        else {
            log->error( "[inputfile] [moment_system] 'quadPoints' not set!" );
            validConfig = false;
        }
        _maxIterations = moment_system->get_as<unsigned>( "maxIterations" ).value_or( 1000 );
        _epsilon       = moment_system->get_as<double>( "epsilon" ).value_or( 5e-5 );
        _hasSource     = false;
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
int Settings::GetWriteFrequency() const { return _writeFrequency; }

// mesh
unsigned Settings::GetMeshDimension() const { return _meshDimension; }
unsigned Settings::GetNumCells() const { return _numCells; }
void Settings::SetNumCells( unsigned n ) {
    _numCells = n;
    _PEforCell.clear();
    _cellIndexPE.clear();
    bool devideChunks = true;
    // int nChunks       = 1;
    if( !devideChunks ) {    // TODO: Fix like quadrature
        int nXPE = int( ( int( _numCells ) - 1 ) / _npes ) + 1;
        if( _mype == _npes - 1 ) {
            nXPE = int( _numCells ) - _mype * nXPE;
            if( nXPE < 0 ) {
                nXPE = 0;
            }
        }
        _nXPE           = unsigned( nXPE );
        unsigned jStart = static_cast<unsigned>( _mype * ( ( static_cast<int>( _numCells ) - 1 ) / _npes + 1.0 ) );
        for( unsigned j = 0; j < _nXPE; ++j ) _cellIndexPE.push_back( jStart + j );
        for( unsigned j = 0; j < _numCells; ++j ) _PEforCell.push_back( int( std::floor( j / _nXPE ) ) );
    }
    else {
        for( unsigned j = 0; j < _numCells; ++j ) {
            _PEforCell.push_back( int( j ) % ( _npes ) );
            if( _PEforCell[j] == _mype ) _cellIndexPE.push_back( j );
        }
        _nXPE = unsigned( _cellIndexPE.size() );
    }
}
std::string Settings::GetOutputFile() const { return _outputFile; }
bool Settings::HasICFile() const { return !_icFile.empty(); }
std::string Settings::GetICFile() const { return _icFile; }
bool Settings::HasRestartFile() const { return !_restartFile.empty(); }
std::string Settings::GetRestartFile() const { return _restartFile; }
bool Settings::HasReferenceFile() const { return !_referenceFile.empty(); }
std::string Settings::GetReferenceFile() const { return _referenceFile; }
bool Settings::LoadLambda() const { return _loadLambda; }

// problem
TimesteppingType Settings::GetTimesteppingType() const { return _timesteppingType; }
double Settings::GetCFL() const { return _CFL; }
double Settings::GetTEnd() const { return _tEnd; }
double Settings::GetMinResidual() const { return _minResidual; }
double Settings::GetDT() const { return _dt; }
void Settings::SetDT( double dt ) { _dt = dt; }
unsigned Settings::GetNDimXi() const { return _numDimXi; }
double Settings::GetGamma() const { return _gamma; }
void Settings::SetGamma( double gamma ) { _gamma = gamma; }
DistributionType Settings::GetDistributionType( unsigned l ) const { return _distributionType[l]; }
std::vector<double> Settings::GetSigma() const { return _sigma; }
double Settings::GetSigma( unsigned l ) const { return _sigma[l]; }
void Settings::SetExactSolution( bool hasExactSolution ) { _hasExactSolution = hasExactSolution; }
bool Settings::HasExactSolution() const { return _hasExactSolution; }
bool Settings::HasSource() const { return _hasSource; }
void Settings::SetSource( bool hasSource ) { _hasSource = hasSource; }

// moment_system
ClosureType Settings::GetClosureType() const { return _closureType; }
void Settings::SetClosureType( ClosureType cType ) { _closureType = cType; }
unsigned Settings::GetNMoments() const { return _nMoments; }
unsigned Settings::GetNQuadPoints() const { return _nQuadPoints; }
void Settings::SetNQuadPoints( unsigned nqNew ) { _nQuadPoints = nqNew; }
unsigned Settings::GetNQTotal() const { return _nQTotal; }
bool Settings::UsesMaxDegree() const { return _useMaxDegree; }
LimiterType Settings::GetLimiterType() const { return _limiterType; }
unsigned Settings::GetMaxIterations() const { return _maxIterations; }
void Settings::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }
double Settings::GetEpsilon() const { return _epsilon; }
unsigned Settings::GetNTotal() const { return _nTotal; }
VectorU Settings::GetNTotalRefinementLevel() const { return _nTotalRefinementLevel; }
std::vector<std::vector<unsigned>> Settings::GetPolyIndices() const { return _polyIndices; }
unsigned Settings::GetNRefinementLevels() const { return _nRefinementLevels; }
unsigned Settings::GetNTotalforRefLevel( unsigned level ) const { return _nTotalRefinementLevel[level]; }
unsigned Settings::GetPolyDegreeforRefLevel( unsigned level ) const { return _refinementLevel[level]; }
GridType Settings::GetGridType() const { return _gridType; }
void Settings::SetNQTotal( unsigned nqTotalNew ) {
    _nQTotal = nqTotalNew;
    // if number of quadrature points is updated, the MPI bounds need to be updated as well
    _kEnd   = unsigned( std::floor( ( double( _mype ) + 1.0 ) * double( _nQTotal / double( _npes ) ) ) );
    _kStart = unsigned( std::ceil( double( _mype ) * ( double( _nQTotal ) / double( _npes ) ) ) );
    if( unsigned( std::ceil( ( double( _mype ) + 1.0 ) * ( double( _nQTotal ) / double( _npes ) ) ) ) == _kEnd ) _kEnd = _kEnd - 1;
    _nQPE = _kEnd - _kStart + 1;
}

// plot
unsigned Settings::GetPlotStepInterval() const { return _plotStepInterval; }
double Settings::GetPlotTimeInterval() const { return _plotTimeInterval; }

// MPI
int Settings::GetMyPE() const { return _mype; }
int Settings::GetNPEs() const { return _npes; }
unsigned Settings::GetKStart() const { return _kStart; }
unsigned Settings::GetKEnd() const { return _kEnd; }
unsigned Settings::GetNqPE() const { return _nQPE; }
unsigned Settings::GetNxPE() const { return _nXPE; }
std::vector<unsigned> Settings::GetCellIndexPE() const { return _cellIndexPE; }
std::vector<int> Settings::GetPEforCell() const { return _PEforCell; }
