#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <string>
#include <iostream>

//#include "mesh.h"
//#include "timediscretization.h"

class Problem
{
private:
    //Mesh* _mesh;
    //TimeDiscretizaiton* _timeDiscretization;
    double _tEnd;
    double _CFL;
    std::string _limiter;
    virtual void solve() = 0;

    // I/O
    std::string _inputFile;
    std::string _outputFolder;
    void parse();
    virtual void print() const = 0;
    virtual void plot() const = 0;
    virtual void writeToFile(std::string filename, int filetype) const = 0;

    Problem(){}

public:
    Problem(std::string inputFile);
    static Problem* create(std::string inputFile);
    ~Problem();
};

#endif // PROBLEM_H
