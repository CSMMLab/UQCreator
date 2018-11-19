#ifndef SETTINGS_H
#define SETTINGS_H

#include <experimental/filesystem>
#include <filesystem>
#include <iostream>
#include <spdlog/spdlog.h>

#include "cpptoml.h"

enum ProblemType { P_BURGERS_1D, P_EULER_1D, P_EULER_2D };
enum ClosureType { C_BOUNDEDBARRIER, C_STOCHASTICGALERKIN, C_EULER_1D, C_EULER_2D };
enum LimiterType { L_MINMOD, L_NONE };
enum TimesteppingType { T_EXPLICITEULER };
enum DistributionType { D_LEGENDRE, D_HERMITE };

class Settings
{
  private:
    // paths
    std::filesystem::path _cwd;
    std::filesystem::path _inputDir;
    std::filesystem::path _inputFile;
    std::filesystem::path _outputDir;
    std::filesystem::path _outputFile;
    std::filesystem::path _continueFile;

    // requied settings
    unsigned _meshDimension;
    unsigned _nQuadPoints;
    unsigned _nMoments;
    unsigned _maxIterations;
    unsigned _nStates;
    unsigned _numCells;
    double _epsilon;
    double _CFL;
    double _tEnd;
    unsigned _plotStepInterval;
    double _plotTimeInterval;
    LimiterType _limiterType;
    ClosureType _closureType;
    ProblemType _problemType;
    TimesteppingType _timesteppingType;
    DistributionType _distributionType;

    // mesh dependent settings

    // problem specific settings
    double _gamma;

    Settings() {}

  public:
    Settings( std::string inputFile );
    ~Settings();

    ProblemType GetProblemType() const;
    unsigned GetNStates() const;
    void SetNStates( unsigned n );
    std::string GetInputFile() const;
    std::string GetInputDir() const;
    std::string GetOutputDir() const;

    // mesh
    unsigned GetMeshDimension() const;
    unsigned GetNumCells() const;
    void SetNumCells( unsigned n );
    std::string GetOutputFile() const;
    bool HasContinueFile() const;
    std::string GetContinueFile() const;

    // problem
    TimesteppingType GetTimesteppingType() const;
    DistributionType GetDistributionType() const;
    double GetCFL() const;
    double GetTEnd() const;
    double GetGamma() const;
    void SetGamma( double gamma );

    // moment_system
    ClosureType GetClosureType() const;
    unsigned GetNMoments() const;
    unsigned GetNQuadPoints() const;
    LimiterType GetLimiterType() const;
    unsigned GetMaxIterations() const;
    double GetEpsilon() const;

    // plot
    bool GetPlotEnabled() const;
    unsigned GetPlotStepInterval() const;
    double GetPlotTimeInterval() const;
};

#endif    // SETTINGS_H
