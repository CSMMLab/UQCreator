#ifndef PROBLEM_H
#define PROBLEM_H

#include <cpptoml.h>
#include <experimental/filesystem>
#include <spdlog/spdlog.h>
#include <string>

#include "mesh.h"
#include "settings.h"
#include "typedefs.h"

class Problem
{
  protected:
    Settings* _settings;
    unsigned _nStates;
    std::shared_ptr<spdlog::logger> _log;
    Mesh* _mesh;

    // variance vector
    Vector _sigma;

    Problem() {}

  public:
    Problem( Settings* settings );
    static Problem* Create( Settings* settings );
    virtual ~Problem();
    virtual void Solve() {}
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) = 0;
    virtual double ComputeDt( const Tensor& u, double dx, unsigned level ) const;
    virtual Vector IC( const Vector& x, const Vector& xi )     = 0;
    virtual Vector LoadIC( const Vector& x, const Vector& xi ) = 0;
    virtual Matrix ExactSolution( double t, const Matrix& x, const Vector& xi ) const;
    virtual Matrix Source( const Matrix& uQ ) const;
    virtual Tensor Source( const Tensor& uQ, const Vector& x, double t, unsigned level ) const;
    virtual void SourceImplicit( Matrix& uQNew, const Matrix& uQTilde, const Matrix& uQ, const Vector& x, double t, unsigned level ) const;
    virtual Matrix BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const;
    virtual void SetMesh( Mesh* mesh ) { _mesh = mesh; }
};

#endif    // PROBLEM_H