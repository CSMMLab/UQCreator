#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/lapack/gesv.h>
#include <blaze/math/lapack/posv.h>
#include <vector>

#include "legendre.h"
#include "settings.h"
#include "typedefs.h"

class Closure
{
  private:
    int* _perm;
    Closure() = delete;

  protected:
    Settings* _settings;
    Polynomial* _basis;
    Polynomial* _quad;
    std::vector<Vector> _phi;    // stores basis functions evaluated at quadrature points
    Matrix _phiTilde;            // stores scaled basis functions evaluated at quadrature points
    Matrix _phiTildeTrans;       // stores scaled basis functions evaluated at quadrature points
    Matrix _phiTildeW;           // stores scaled basis functions evaluated at quadrature points times weight
    std::vector<Vector> _phiTildeVec;
    std::vector<Matrix> _hPartial;    // stores partial matrices for Hessian computation
    double _alpha;
    double _uMinus, _uPlus;    // IPM bounds for scalar problems
    unsigned _nMoments;
    unsigned _nQuadPoints;
    unsigned _nStates;
    void Hessian( Matrix& H, const Matrix& lambda );
    void Gradient( Vector& g, const Matrix& lambda, const Matrix& u );

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
    virtual Matrix SolveClosure( const Matrix& uMatrix, Matrix& lambda );
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    Vector EvaluateLambda( const Matrix& lambda, unsigned k );
    Matrix EvaluateLambda( const Matrix& lambda ) const;
    Vector EvaluateLambda( const Matrix& lambda, const Vector& xi, unsigned k );
    Matrix EvaluateLambda( const Matrix& lambda, const Vector& xi );
    /**
     * calculate solution for kinetic entropy with given entropic variable
     */
    virtual void U( Vector& out, const Vector& Lambda ) = 0;
    virtual Matrix U( const Matrix& Lambda )            = 0;
    /**
     * calculate derivative of solution for kinetic entropy with given entropic variable
     */
    virtual void DU( Matrix& y, const Vector& Lambda ) = 0;

    const std::vector<Vector>& GetPhi() { return _phi; }
    const Vector& GetPhiTilde( int k ) { return _phiTildeVec[k]; }
    const Matrix& GetPhiTilde() { return _phiTilde; }
    const Matrix& GetPhiTildeW() { return _phiTildeW; }
    /**
     * Transform matrix to vector
     * @return
     */
    Vector MakeVector( const Matrix& mat ) const;
    /**
     * Transform vector to matrix
     * @return
     */
    Matrix MakeMatrix( const Vector& vec ) const;

    double GetUPlus() const { return _uPlus; }
    double GetUMinus() const { return _uMinus; }
    double CalcNorm( Vector& test );
};

#endif    // CLOSURE_H
