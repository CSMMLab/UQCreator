#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <string>

#include "mesh.h"
//#include "timediscretization.h"

#define QUAD_TYPE_LEGENDRE 201
#define QUAD_TYPE_HERMITE 202

class Problem
{
protected:
    Mesh* _mesh;
//    TimeDiscretization* _timeDiscretization;
    int _quadType;
    int _nQuadPoints;
    int _nMoments;
    int _maxIterations;
    double _epsilon;
    double _CFL;
    double _tEnd;
    std::string _limiter;

    // I/O
    std::string _inputFile;
    std::string _outputFolder;
    virtual void Print() const = 0;
    virtual void WriteToFile(std::string filename, int filetype) const = 0;

    Problem(){}

public:
    Problem(std::string inputFile);
    static Problem* Create(std::string inputFile);
    virtual ~Problem();
    virtual void Solve(){}
    virtual double G(double u, double v) = 0;
    virtual void Plot(blaze::DynamicVector<double>& x, blaze::DynamicVector<double>& u) const = 0;

    int GetQuadType();
    int GetNQuadPoints();
    int GetNMoments();
    int GetMaxIterations();
    double GetEpsilon();
    double GetCFL();
    double GetTEnd();
    Mesh* GetMesh();
};

#endif // PROBLEM_H
