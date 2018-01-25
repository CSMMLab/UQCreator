#ifndef CLOSURE_H
#define CLOSURE_H

#include "legendre.h"
#include "problem.h"
#include <blaze/math/LAPACK.h>
#include <vector>

class Closure
{
  protected:
    Problem* _problem;
    Polynomial* _basis;
    Polynomial* _quad;
    std::vector<blaze::DynamicVector<double>> _phi;    // stores basis functions evaluated at quadrature points
    blaze::DynamicMatrix<double> _phiTilde;            // stores scaled basis functions evaluated at quadrature points
    blaze::DynamicMatrix<double> _phiTildeTrans;       // stores scaled basis functions evaluated at quadrature points
    blaze::DynamicMatrix<double> _phiTildeW;           // stores scaled basis functions evaluated at quadrature points times weight
    std::vector<blaze::DynamicVector<double>> _phiTildeVec;
    std::vector<blaze::DynamicMatrix<double>> _hPartial;    // stores partial matrices for Hessian computation
    double _uMinus, _uPlus;                                 // IPM bounds for scalar problems
    int _nMoments;
    int _nQuadPoints;
    int _nStates;
    void Hessian( blaze::DynamicMatrix<double>& H, const blaze::DynamicMatrix<double>& lambda );
    void Gradient( blaze::DynamicVector<double>& g, const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicMatrix<double>& u );

  public:
    /**
     * constructor of class Closure
     * @param pointer to problem class
     */
    Closure( Problem* problem );
    static Closure* Create( Problem* problem );
    /**
     * calculate dual vector fulfilling the moment constraint
     * @param moment vector
     * @param initial guess for dual vector
     * @return correct dual vector
     */
    blaze::DynamicMatrix<double> SolveClosure( const blaze::DynamicMatrix<double>& uMatrix, blaze::DynamicMatrix<double>& lambda );
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    blaze::DynamicVector<double> EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, int k );
    blaze::DynamicMatrix<double> EvaluateLambda( const blaze::DynamicMatrix<double>& lambda ) const;
    blaze::DynamicVector<double> EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicVector<double>& xi, int k );
    blaze::DynamicMatrix<double> EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicVector<double>& xi );
    /**
     * calculate solution for kinetic entropy with given entropic variable
     */
    virtual void U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda ) = 0;
    virtual blaze::DynamicMatrix<double> U( const blaze::DynamicMatrix<double>& Lambda )            = 0;
    /**
     * calculate derivative of solution for kinetic entropy with given entropic variable
     */
    virtual void DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda ) = 0;

    const std::vector<blaze::DynamicVector<double>>& GetPhi() { return _phi; }
    const blaze::DynamicVector<double>& GetPhiTilde( int k ) { return _phiTildeVec[k]; }
    const blaze::DynamicMatrix<double>& GetPhiTilde() { return _phiTilde; }
    const blaze::DynamicMatrix<double>& GetPhiTildeW() { return _phiTildeW; }
    /**
     * Transform matrix to vector
     * @return
     */
    blaze::DynamicVector<double> MakeVector( const blaze::DynamicMatrix<double>& mat ) const;
    /**
     * Transform vector to matrix
     * @return
     */
    blaze::DynamicMatrix<double> MakeMatrix( const blaze::DynamicVector<double>& vec ) const;

    double GetUPlus() const { return _uPlus; }
    double GetUMinus() const { return _uMinus; }
    double CalcNorm( blaze::DynamicVector<double>& test );
};

#endif    // CLOSURE_H
