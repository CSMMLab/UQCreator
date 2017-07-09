#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <string>
#include <iostream>

#include "mesh.h"
//#include "timediscretization.h"

#define QUAD_TYPE_LEGENDRE 201
#define QUAD_TYPE_HERMITE 202

class Problem
{
protected:
    Mesh* _mesh;
    //TimeDiscretizaiton* _timeDiscretization;
    int _quadType;
    int _nQuadPoints;
    virtual void solve() = 0;

    // I/O
    std::string _inputFile;
    std::string _outputFolder;
    virtual void print() const = 0;
    virtual void plot() const = 0;
    virtual void writeToFile(std::string filename, int filetype) const = 0;

    Problem(){}

public:
    Problem(std::string inputFile);
    static Problem* create(std::string inputFile);
    virtual ~Problem();

    int getQuadType();
    int getNQuadPoints();
};

#endif // PROBLEM_H
