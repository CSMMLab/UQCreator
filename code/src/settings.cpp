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
            else if( problemTypeString->compare( "PNEquations" ) == 0 ) {
                if( _meshDimension == 1 )
                    _problemType = ProblemType::P_PNEQUATIONS_1D;
                else if( _meshDimension == 2 )
                    _problemType = ProblemType::P_PNEQUATIONS_2D;
            }
            else if( problemTypeString->compare( "M1Equations" ) == 0 ) {
                if( _meshDimension == 1 )
                    _problemType = ProblemType::P_M1EQUATIONS_1D;
                else
                    std::cerr << "M1 for 2D not implemented" << std::endl;
            }
            else if( problemTypeString->compare( "RadiationHydrodynamics" ) == 0 ) {
                if( _meshDimension == 1 )
                    _problemType = ProblemType::P_RADIATIONHYDRO_1D;
                else if( _meshDimension == 2 )
                    _problemType = ProblemType::P_RADIATIONHYDRO_2D;
            }
            else if( problemTypeString->compare( "ThermalRadiativeTransfer" ) == 0 ) {
                _problemType = ProblemType::P_THERMALRAD_1D;
            }
            else if( problemTypeString->compare( "ThermalPN" ) == 0 ) {
                _problemType = ProblemType::P_THERMALPN_1D;
            }
            else if( problemTypeString->compare( "NavierStokes" ) == 0 ) {
                _problemType = ProblemType::P_NAVIERSTOKES_1D;
            }
            else {
                log->error( "[inputfile] [general] 'problem' is invalid!\nPlease set one of the following types: Burgers1D, Euler1D, "
                            "Euler2D,ShallowWater1D,ShallowWater2D,PNEquations2D,RadiationHydrodynamics,M1Equations,NavierStokes" );
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
        auto outputFile = mesh->get_as<std::string>( "outputFile" );
        if( outputFile ) {
            _outputFile = _inputDir.string() + "/" + _outputDir.string() + "/" + *outputFile;
        }
        else {
            log->error( "[inputfile] [mesh] 'outputFile' not set!" );
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
            _writeFrequency = general->get_as<int>( "writeFrequency" ).value_or( -1 );
            if( _writeFrequency == -1 )
                _writeInTime = false;
            else
                _writeInTime = true;
            if( refFile ) {
                _referenceFile = _inputDir.string() + "/" + *refFile;
            }
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
            else if( closureTypeString->compare( "RadiationHydrodynamics" ) == 0 ) {
                _closureType = ClosureType::C_RADHYDRO;
            }
            else if( closureTypeString->compare( "ThermalRadiation" ) == 0 ) {
                _closureType = ClosureType::C_THERMALRAD_1D;
            }
            else if( closureTypeString->compare( "M1" ) == 0 ) {
                _closureType = ClosureType::C_M1_1D;
            }
            else if( closureTypeString->compare( "HyperbolicityLimiter" ) == 0 ) {
                _closureType = ClosureType::C_HYPLIM;
            }
            else {
                log->error( "[inputfile] [moment_system] 'closure' is invalid!\nPlease set one of the following types: BoundedBarrier, LogSin, "
                            "StochasticGalerkin, Euler, Euler2D,L2Filter,LassoFilter,RadiationHydrodynamics,M1,HyperbolicityLimiter" );
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

            _polyIndices.resize( 0 );
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
                    }
                }
                _nTotalRefinementLevel[level] = static_cast<unsigned>( _polyIndices.size() );
                if( UsesMaxDegree() ) break;
                previousDegree = int( GetPolyDegreeforRefLevel( level ) );
            }
            _nTotal = unsigned( _polyIndices.size() );
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
            assert( nQArray->size() == _nRefinementLevels );
            _quadLevel = VectorU( _nRefinementLevels );
            for( unsigned i = 0; i < _nRefinementLevels; ++i ) _quadLevel[i] = unsigned( ( *nQArray )[i] );

            // computing Nq Total
            if( quadratureType->at( 0 ).compare( "sparseGrid" ) == 0 ) {
                _gridType = G_SPARSEGRID;
                _nQTotal  = unsigned( std::pow( 2, _nQuadPoints ) ) + 1u;
            }
            else if( quadratureType->at( 0 ).compare( "tensorizedGrid" ) == 0 ) {    // tensorizedGrid
                _gridType = G_TENSORIZEDGRID;
                _nQTotal  = unsigned( std::pow( _nQuadPoints, _numDimXi ) );
            }
            else if( quadratureType->at( 0 ).compare( "tensorizedCCGrid" ) == 0 ) {    // tensorizedCC
                _gridType = G_TENSORIZEDCC;
                if( _nQuadPoints == 0 )
                    _nQTotal = 1;
                else
                    _nQTotal = unsigned( std::pow( static_cast<unsigned>( std::pow( 2, _nQuadPoints ) + 1 ), _numDimXi ) );
            }
            else {
                log->error( "[inputfile] [moment_system] 'quadType' not defined!" );
                validConfig = false;
            }
        }
        else {
            log->error( "[inputfile] [moment_system] 'quadPoints' not set!" );
            validConfig = false;
        }

        // read convergence retardation
        auto cRetardation = moment_system->get_array_of<cpptoml::array>( "cRetardation" );
        if( cRetardation ) {
            auto retardationSteps    = ( *cRetardation )[0]->get_array_of<int64_t>();
            auto residualRetardation = ( *cRetardation )[1]->get_array_of<double>();
            assert( retardationSteps->size() == residualRetardation->size() );
            _nRetardationLevels  = unsigned( retardationSteps->size() ) + 1;
            _retardationSteps    = VectorU( _nRetardationLevels, false );
            _residualRetardation = Vector( _nRetardationLevels, false );

            for( unsigned i = 0; i < _nRetardationLevels - 1; ++i ) {
                _retardationSteps[i]    = unsigned( ( *retardationSteps )[i] );
                _residualRetardation[i] = ( *residualRetardation )[i];
            }
            _retardationSteps[_nRetardationLevels - 1]    = _nRefinementLevels;
            _residualRetardation[_nRetardationLevels - 1] = _minResidual;
        }
        else {
            _retardationSteps    = VectorU( 1, _nRefinementLevels );
            _residualRetardation = Vector( 1, _minResidual );
            _nRetardationLevels  = 1;
        }
        if( _nRefinementLevels > 1 ) {
            auto refinementThresholds = moment_system->get_array_of<double_t>( "refinementThresholds" );
            // assert( refinementThresholds->size() == 2 );
            if( refinementThresholds->size() != 2 ) {
                std::cerr << "[settings]: No Refinement Barrier Specified." << std::endl;
                exit( EXIT_FAILURE );
            }
            _refinementThreshold = ( *refinementThresholds )[0];
            _coarsenThreshold    = ( *refinementThresholds )[1];
        }

        _regularizationStrength = moment_system->get_as<double>( "regularizationStrength" ).value_or( -1.0 );
        if( _regularizationStrength < 0 )
            _regularization = false;
        else
            _regularization = true;

        _filterStrength = moment_system->get_as<double>( "filterStrength" ).value_or( -1.0 );

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
unsigned Settings::GetMaxIterations() const { return _maxIterations; }
void Settings::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }
double Settings::GetEpsilon() const { return _epsilon; }
unsigned Settings::GetNTotal() const { return _nTotal; }
VectorU Settings::GetNTotalRefinementLevel() const { return _nTotalRefinementLevel; }
std::vector<std::vector<unsigned>> Settings::GetPolyIndices() const { return _polyIndices; }
unsigned Settings::GetNRefinementLevels() const { return _nRefinementLevels; }
unsigned Settings::GetNRefinementLevels( unsigned retardation ) const { return _retardationSteps[retardation]; }
double Settings::GetRefinementThreshold() const { return _refinementThreshold; }
void Settings::SetRefinementThreshold( double ref ) { _refinementThreshold = ref; }
double Settings::GetCoarsenThreshold() const { return _coarsenThreshold; }
void Settings::SetCoarsenThreshold( double coa ) { _coarsenThreshold = coa; }
double Settings::GetResidualRetardation( unsigned retardation ) const { return _residualRetardation[retardation]; }
unsigned Settings::GetNRetardationLevels() const { return _nRetardationLevels; }
unsigned Settings::GetNTotalforRefLevel( unsigned level ) const { return _nTotalRefinementLevel[level]; }
unsigned Settings::GetPolyDegreeforRefLevel( unsigned level ) const { return _refinementLevel[level]; }
GridType Settings::GetGridType() const { return _gridType; }
VectorU Settings::GetQuadLevel() const { return _quadLevel; }
std::vector<unsigned> Settings::GetIndicesQforRef( unsigned level ) const { return _kIndicesAtRef[level]; }
bool Settings::HasRegularization() const { return _regularization; }

// Set Total number of Quadrature points at each refinement level
void Settings::SetNQTotalForRef( const VectorU& nQTotalForRef ) {
    _nQTotalForRef = nQTotalForRef;
    _nQPEAtRef     = VectorU( _nRefinementLevels );
    unsigned kEnd, kStart;
    _kIndicesAtRef.resize( _nRefinementLevels );
    unsigned nQTotalForRefOld = 0;
    for( unsigned l = 0; l < _nRefinementLevels; ++l ) {
        unsigned numberNewPoints = _nQTotalForRef[l] - nQTotalForRefOld;
        // compute end and start point for each refinement level for standard distribution strategy
        kEnd   = unsigned( std::floor( ( double( _mype ) + 1.0 ) * double( double( numberNewPoints ) / double( _npes ) ) ) );
        kStart = unsigned( std::ceil( double( _mype ) * ( double( numberNewPoints ) / double( _npes ) ) ) );
        if( unsigned( std::ceil( ( double( _mype ) + 1.0 ) * ( double( numberNewPoints ) / double( _npes ) ) ) ) == kEnd ) kEnd = kEnd - 1;
        kEnd += nQTotalForRefOld;
        kStart += nQTotalForRefOld;

        // clear vector
        _kIndicesAtRef[l].clear();
        _kIndicesAtRef[l].resize( 0 );

        // save old indices on current refinement level
        if( l > 0 )
            for( unsigned i = 0; i < _kIndicesAtRef[l - 1].size(); ++i ) _kIndicesAtRef[l].push_back( _kIndicesAtRef[l - 1][i] );

        // save new indices on current refinement level
        for( unsigned k = kStart; k <= kEnd; ++k ) {
            _kIndicesAtRef[l].push_back( k );
        }

        _nQPEAtRef[l]    = unsigned( _kIndicesAtRef[l].size() );
        nQTotalForRefOld = _nQTotalForRef[l];
    }
    _nQPE    = _nQPEAtRef[_nRefinementLevels - 1];
    _nQTotal = _nQTotalForRef[_nRefinementLevels - 1];
}
unsigned Settings::GetNQTotalForRef( unsigned level ) const { return _nQTotalForRef[level]; }
VectorU Settings::GetNQTotalForRef() const { return _nQTotalForRef; }
unsigned Settings::GetNqPEAtRef( unsigned level ) const { return _nQPEAtRef[level]; }
double Settings::GetFilterStrength() const { return _filterStrength; }
double Settings::GetRegularizationStrength() const { return _regularizationStrength; }

// plot
unsigned Settings::GetPlotStepInterval() const { return _plotStepInterval; }
double Settings::GetPlotTimeInterval() const { return _plotTimeInterval; }

// MPI
int Settings::GetMyPE() const { return _mype; }
int Settings::GetNPEs() const { return _npes; }
unsigned Settings::GetNqPE() const { return _nQPE; }
unsigned Settings::GetNxPE() const { return _nXPE; }
std::vector<unsigned> Settings::GetCellIndexPE() const { return _cellIndexPE; }
std::vector<int> Settings::GetPEforCell() const { return _PEforCell; }
