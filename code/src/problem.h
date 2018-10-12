#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <string>

#include "mesh.h"
#include "typedefs.h"

#define QUAD_TYPE_LEGENDRE 201
#define QUAD_TYPE_HERMITE 202

class Problem
{
  protected:
    Mesh* _mesh;
    int _quadType;
    unsigned _nQuadPoints;
    unsigned _nMoments;
    unsigned _maxIterations;
    unsigned _nStates;
    double _epsilon;
    double _CFL;
    double _tEnd;
    std::string _limiter;
    std::string _closureType;
    std::string _problemType;

    // I/O
    std::string _inputFile;
    std::string _outputDir;

    Problem() {}

  public:
    Problem( std::string inputFile );
    static Problem* Create( std::string inputFile );
    virtual ~Problem();
    virtual void Solve() {}
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) = 0;
    virtual double ExactSolution( double t, double x, double xi )                              = 0;
    virtual double GetGamma() const { return -1.0; }

    int GetQuadType() const;
    unsigned GetNQuadPoints() const;
    unsigned GetNMoments() const;
    unsigned GetMaxIterations() const;
    unsigned GetNStates() const;
    double GetEpsilon() const;
    double GetCFL() const;
    double GetTEnd() const;
    Mesh* GetMesh() const;
    std::string GetInputFile() const;
    std::string GetLimiter() const;
    std::string GetOutputDir() const;
    std::string GetClosureType() const;
    std::string GetProblemType() const;
};

#endif    // PROBLEM_H
