#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <string>
#include <iostream>

class Problem
{
private:
    int _dimension;
    int* _discretization;
    double _tEnd;
    double _CFL;
    std::string _limiter;

    std::string _inputFile;
    std::string _outputFolder;

    void parse();

    Problem(){}

public:
    Problem(std::string inputFile);
    static Problem* create(std::string inputFile);
    ~Problem();
};

#endif // PROBLEM_H
