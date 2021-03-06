#ifndef SETTINGS_H
#define SETTINGS_H

#include <experimental/filesystem>
#include <filesystem>
#include <iostream>
#include <spdlog/spdlog.h>

#include "cpptoml.h"

#include "typedefs.h"

enum ProblemType { P_BURGERS_1D, P_EULER_1D, P_EULER_2D, P_SHALLOWWATER_1D, P_SHALLOWWATER_2D, P_PNEQUATIONS_2D };
enum ClosureType {
    C_BOUNDEDBARRIER,
    C_LOGSIN,
    C_STOCHASTICGALERKIN,
    C_EULER_1D,
    C_EULER_2D,
    C_SHALLOWWATER_1D,
    C_SHALLOWWATER_2D,
    C_L2FILTER,
    C_LASSOFILTER,
    C_EXPFILTER
};
enum LimiterType { L_MINMOD, L_NONE };
enum TimesteppingType { T_EXPLICITEULER };
enum DistributionType { D_LEGENDRE, D_HERMITE };
enum GridType { G_SPARSEGRID, G_TENSORIZEDGRID };

class Settings
{
  private:
    // paths
    std::filesystem::path _cwd;
    std::filesystem::path _inputDir;
    std::filesystem::path _inputFile;
    std::filesystem::path _outputDir;
    std::filesystem::path _outputFile;
    std::filesystem::path _icFile;
    std::filesystem::path _restartFile;
    std::filesystem::path _referenceFile;
    bool _loadLambda;

    int _writeFrequency;    // number of time steps until error to reference solution is computed

    // requied settings
    unsigned _meshDimension;
    unsigned _numCells;
    unsigned _nXPE;    // number of spatial cells on PE

    unsigned _nQuadPoints;    // number of quadrature points in one dimension
    unsigned _nQTotal;        // number of quadrature points in all dimensions
    unsigned _nQPE;           // number of total quadrature points on PE
    unsigned _kStart;         // start point in quadrature point array for PE
    unsigned _kEnd;           // end point in quadrature point array for PE
    VectorU _quadLevel;       // quadrature level array

    int _mype;                             // PE number
    int _npes;                             // number of all PEs
    std::vector<unsigned> _cellIndexPE;    // vector of spatial cells for PE
    std::vector<int> _PEforCell;

    unsigned _nMoments;                                 // number of moments in one dimension
    unsigned _nTotal;                                   // number of moments in all dimensions
    VectorU _nTotalRefinementLevel;                     // vector of number of moments in all dimensions for each refinement level
    VectorU _refinementLevel;                           // vector of different refinement levels
    std::vector<std::vector<unsigned>> _polyIndices;    // indices for polynomial basis functions
    unsigned _nRefinementLevels;
    bool _useMaxDegree;    // specifies moment hierarchy

    unsigned _maxIterations;
    unsigned _nStates;

    bool _hasExactSolution;    // indicates if exact solution is implemented in Problem

    unsigned _numDimXi;
    double _epsilon;
    double _CFL;
    double _tEnd;
    double _minResidual;    // residual at which iteration is stopped for steady problems
    double _dt;             // timestepsize only required if no function ComputeDt provided in problem
    unsigned _plotStepInterval;
    double _plotTimeInterval;
    LimiterType _limiterType;
    ClosureType _closureType;
    ProblemType _problemType;
    TimesteppingType _timesteppingType;
    GridType _gridType;
    std::vector<DistributionType> _distributionType;
    std::vector<double> _sigma;

    // mesh dependent settings

    // problem specific settings
    double _gamma;
    bool _hasSource;

    Settings() = delete;
    void Init( std::shared_ptr<cpptoml::table> file, bool restart );

  public:
    Settings( std::string inputFile );
    Settings( const std::istringstream& inputStream );
    ~Settings();

    ProblemType GetProblemType() const;
    unsigned GetNStates() const;
    void SetNStates( unsigned n );
    std::string GetInputFile() const;
    std::string GetInputDir() const;
    std::string GetOutputDir() const;
    int GetWriteFrequency() const;

    // mesh
    unsigned GetMeshDimension() const;
    unsigned GetNumCells() const;
    void SetNumCells( unsigned n );
    std::string GetOutputFile() const;
    bool HasICFile() const;
    std::string GetICFile() const;
    bool HasRestartFile() const;
    std::string GetRestartFile() const;
    bool HasReferenceFile() const;
    std::string GetReferenceFile() const;
    bool LoadLambda() const;

    // problem
    TimesteppingType GetTimesteppingType() const;
    DistributionType GetDistributionType( unsigned l ) const;
    std::vector<double> GetSigma() const;
    double GetSigma( unsigned l ) const;
    double GetCFL() const;
    unsigned GetNDimXi() const;
    double GetTEnd() const;
    double GetMinResidual() const;
    double GetDT() const;
    void SetDT( double dt );
    double GetGamma() const;
    void SetGamma( double gamma );
    void SetExactSolution( bool hasExactSolution );
    bool HasExactSolution() const;
    bool HasSource() const;
    void SetSource( bool hasSource );

    // moment_system
    ClosureType GetClosureType() const;
    void SetClosureType( ClosureType cType );
    unsigned GetNMoments() const;

    unsigned GetNTotal() const;
    VectorU GetNTotalRefinementLevel() const;
    VectorU GetQuadLevel() const;
    std::vector<std::vector<unsigned>> GetPolyIndices() const;
    unsigned GetNRefinementLevels() const;
    unsigned GetNTotalforRefLevel( unsigned level ) const;
    unsigned GetPolyDegreeforRefLevel( unsigned level ) const;
    unsigned GetNQuadPoints() const;
    void SetNQuadPoints( unsigned nqNew );
    void SetNQTotal( unsigned nqTotalNew );
    unsigned GetNQTotal() const;
    bool UsesMaxDegree() const;
    LimiterType GetLimiterType() const;
    unsigned GetMaxIterations() const;
    void SetMaxIterations( unsigned maxIterations );
    double GetEpsilon() const;
    GridType GetGridType() const;

    // plot
    bool GetPlotEnabled() const;
    unsigned GetPlotStepInterval() const;
    double GetPlotTimeInterval() const;

    // MPI
    int GetMyPE() const;
    int GetNPEs() const;
    unsigned GetKStart() const;
    unsigned GetKEnd() const;
    unsigned GetNqPE() const;
    unsigned GetNxPE() const;
    std::vector<unsigned> GetCellIndexPE() const;
    std::vector<int> GetPEforCell() const;
};

#endif    // SETTINGS_H
