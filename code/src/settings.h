#ifndef SETTINGS_H
#define SETTINGS_H

#include <experimental/filesystem>
#include <iostream>

#include "cpptoml.h"

enum ProblemType { P_BURGERS_1D, P_EULER_1D, P_EULER_2D };
enum ClosureType { C_BOUNDEDBARRIER, C_STOCHASTICGALERKIN, C_EULER_1D, C_EULER_2D };
enum LimiterType { L_MINMOD, L_NONE };
enum TimesteppingType { T_EXPLICITEULER, T_HEUN, T_SSPMULTISTEP };

class Settings
{
  private:
    std::string _inputFile;
    std::experimental::filesystem::path _inputDir;
    std::experimental::filesystem::path _outputDir;

    // requied settings
    unsigned _nQuadPoints;
    unsigned _nMoments;
    unsigned _maxIterations;
    unsigned _nStates;
    unsigned _numCells;
    double _epsilon;
    double _CFL;
    double _tEnd;
    LimiterType _limiter;
    ClosureType _closureType;
    ProblemType _problemType;

    // mesh dependent settings

    Settings() {}

  public:
    Settings( std::string inputFile );
    ~Settings();

    int GetQuadType() const;
    unsigned GetNQuadPoints() const;
    unsigned GetNMoments() const;
    unsigned GetMaxIterations() const;
    unsigned GetNStates() const;
    double GetEpsilon() const;
    double GetCFL() const;
    double GetTEnd() const;
    std::string GetInputFile() const;
    std::string GetOutputDir() const;
    LimiterType GetLimiter() const;
    ClosureType GetClosureType() const;
    ProblemType GetProblemType() const;
    unsigned GetNumCells() const;
    void SetNumCells( unsigned n ) const;
};

#endif    // SETTINGS_H
