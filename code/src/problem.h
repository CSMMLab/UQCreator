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

    // variance vector
    Vector _sigma;

    Problem() {}

  public:
    Problem( Settings* settings );
    static Problem* Create( Settings* settings );
    virtual ~Problem();
    virtual void Solve() {}
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) = 0;
    virtual double ComputeDt( const Matrix& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi )     = 0;
    virtual Vector LoadIC( const Vector& x, const Vector& xi ) = 0;
    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
    virtual Matrix Source( const Matrix& uQ ) const;
    virtual Matrix Source( const Matrix& uQ, const Vector& x, double t ) const;
    virtual Matrix BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const;
};

#endif    // PROBLEM_H
