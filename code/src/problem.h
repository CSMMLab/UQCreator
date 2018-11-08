#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <experimental/filesystem>
#include <spdlog/spdlog.h>
#include <string>

#include "settings.h"
#include "typedefs.h"

class Problem
{
  protected:
    Settings* _settings;
    unsigned _nStates;
    std::shared_ptr<spdlog::logger> _log;

    Problem() {}

  public:
    Problem( Settings* settings );
    static Problem* Create( Settings* settings );
    virtual ~Problem();
    virtual void Solve() {}
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) = 0;
    virtual double ComputeDt( Vector& u, double dx ) const                                     = 0;
};

#endif    // PROBLEM_H
