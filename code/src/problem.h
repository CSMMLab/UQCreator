#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <experimental/filesystem>
#include <string>

#include "settings.h"
#include "typedefs.h"

class Problem
{
  protected:
    const Settings* _settings;
    unsigned _nStates;

    Problem() {}

  public:
    Problem( const Settings* settings );
    static Problem* Create( const Settings* settings );
    virtual ~Problem();
    virtual void Solve() {}
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) = 0;
    virtual double ExactSolution( double t, double x, double xi )                              = 0;
    virtual double GetGamma() const { return -1.0; }
};

#endif    // PROBLEM_H
