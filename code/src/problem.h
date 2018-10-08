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

    // I/O
    std::string _inputFile;
    std::string _outputDir;
    virtual void Print()                                                 = 0;
    virtual void WriteToFile( std::string filename, int filetype ) const = 0;

    Problem() {}

  public:
    Problem( std::string inputFile );
    static Problem* Create( std::string inputFile );
    virtual ~Problem();
    virtual void Solve() {}
    virtual double G( double u, double v )                        = 0;
    virtual Matrix G( const Matrix& u, const Matrix& v )          = 0;
    virtual void Plot( Vector& x, Vector& u )                     = 0;
    virtual double ExactSolution( double t, double x, double xi ) = 0;
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
};

#endif    // PROBLEM_H
