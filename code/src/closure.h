#ifndef CLOSURE_H
#define CLOSURE_H

#include <spdlog/spdlog.h>
#include <vector>

#include "legendre.h"
#include "settings.h"
#include "typedefs.h"

#include "quadraturegrid.h"

class Closure
{
  private:
    Closure() = delete;

  protected:
    Settings* _settings;
    std::vector<Polynomial*> _basis;
    std::vector<Polynomial*> _quad;
    QuadratureGrid* _quadGrid;
    std::vector<Vector> _xiGrid;
    std::vector<Vector> _wGrid;
    std::vector<Vector> _quadNodes;
    VectorU _nQTotalForRef;
    VectorU _nTotalForRef;
    Matrix _phiTilde;         // stores scaled basis functions evaluated at quadrature points
    Matrix _phiTildeTrans;    // stores scaled basis functions evaluated at quadrature points
    Matrix _phiTildeF;        // stores scaled basis functions evaluated at quadrature points times pdf
    Matrix _phiTildeWf;       // stores scaled basis functions evaluated at quadrature points times weight and pdf
    std::vector<Vector> _phiTildeVec;
    MatVec _hPartial;    // stores partial matrices for Hessian computation
    double _alpha;       // step size for Newton
    unsigned _nMoments;
    unsigned _nQuadPoints;
    unsigned _nStates;
    unsigned _numDimXi;
    unsigned _nQTotal;
    unsigned _nTotal;
    unsigned _maxIterations;
    void Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel );
    void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel );
    std::shared_ptr<spdlog::logger> _log;
    Matrix _dUdLambda;    // preallocated memory dor computation of Hessian
  public:
    /**
     * constructor of class Closure
     * @param pointer to problem class
     */
    Closure( Settings* settings );
    virtual ~Closure();
    static Closure* Create( Settings* settings );
    /**
     * calculate dual vector fulfilling the moment constraint
     * @param moment vector
     * @param initial guess for dual vector
     * @return correct dual vector
     */
    // virtual void SolveClosure( Matrix& lambda, const Matrix& u );
    virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
    virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    Vector EvaluateLambda( const Matrix& lambda, unsigned k, unsigned nTotal );
    Matrix EvaluateLambda( const Matrix& lambda ) const;
    void EvaluateLambda( Matrix& out, const Matrix& lambda ) const;
    Vector EvaluateLambda( const Matrix& lambda, const Vector& xi, unsigned k );
    Matrix EvaluateLambda( const Matrix& lambda, const Vector& xi );
    Matrix EvaluateLambdaOnPE( const Matrix& lambda, unsigned levelOld, unsigned levelNew ) const;
    /**
     * calculate solution for kinetic entropy with given entropic variable
     */
    virtual void U( Vector& out, const Vector& Lambda ) = 0;
    virtual void U( Matrix& out, const Matrix& Lambda ) = 0;
    virtual Matrix U( const Matrix& Lambda )            = 0;
    /**
     * calculate derivative of solution for kinetic entropy with given entropic variable
     */
    virtual void DU( Matrix& y, const Vector& Lambda ) = 0;

    virtual void DS( Vector& ds, const Vector& u ) const;

    const Vector& GetPhiTilde( int k ) const { return _phiTildeVec[k]; }
    const Matrix& GetPhiTilde() const { return _phiTilde; }
    const Matrix& GetPhiTildeWf() const { return _phiTildeWf; }
    const Matrix GetPhiTildeWfAtRef( unsigned level ) const;
    /**
     * Add matrix A and vector b and save result in a matrix
     */
    void AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const;
    /**
     * Add matrix A and vector b and save result in b
     */
    void SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const;

    double CalcNorm( Vector& test, unsigned nTotal ) const;

    /**
     * reset step size for Newton
     * @param new step size
     */
    void SetAlpha( double alpha );

    /**
     * reset max iteration number for Newton
     * @param new max iteration number
     */
    void SetMaxIterations( unsigned maxIterations );

    /**
     * get max iteration number for Newton
     * @return max iteration number
     */
    unsigned GetMaxIterations() const;

    std::vector<Polynomial*> GetBasis();
    std::vector<Polynomial*> GetQuadrature();
    QuadratureGrid* GetQuadratureGrid();
};

#endif    // CLOSURE_H
