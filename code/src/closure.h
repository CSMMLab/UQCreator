#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/LAPACK.h>
#include <vector>
#include "problem.h"
#include "legendre.h"

class Closure
{
private:
    Problem* _problem;
    Polynomial* _basis;
    Polynomial* _quad;
    std::vector<blaze::DynamicVector<double>> _phi; // stores basis functions evaluated at quadrature points
    std::vector<blaze::DynamicVector<double>> _phiTilde; // stores scaled basis functions evaluated at quadrature points
    std::vector<blaze::DynamicMatrix<double>> _hPartial; // stores partial matrices for Hessian computation
    double _uMinus, _uPlus; // IPM bounds for scalar problems
    int _nMoments;
    int _nQuadPoints;
    int _nStates;
    blaze::DynamicMatrix<double> Hessian(const blaze::DynamicMatrix<double> &lambda);
public:
    blaze::DynamicVector<double> Gradient(const blaze::DynamicMatrix<double> &lambda, const blaze::DynamicMatrix<double> &u);
public:
    /**
     * constructor of class Closure
     * @param pointer to problem class
     */
    Closure(Problem *problem);
    /**
     * calculate dual vector fulfilling the moment constraint
     * @param moment vector
     * @param initial guess for dual vector
     * @return correct dual vector
     */
    blaze::DynamicMatrix<double> SolveClosure(const blaze::DynamicMatrix<double> &uMatrix, blaze::DynamicMatrix<double> &lambda);
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    blaze::DynamicVector<double> EvaluateLambda(const blaze::DynamicMatrix<double>& lambda, int k);
    blaze::DynamicMatrix<double> EvaluateLambda(const blaze::DynamicMatrix<double>& lambda)const;
    blaze::DynamicVector<double> EvaluateLambda(const blaze::DynamicMatrix<double> &lambda, const blaze::DynamicVector<double>& xi, int k);
    blaze::DynamicMatrix<double> EvaluateLambda(const blaze::DynamicMatrix<double> &lambda, const blaze::DynamicVector<double>& xi);
    /**
     * calculate solution for kinetic entropy with given entropic variable
     * @param entropic variable
     * @return solution
     */
    blaze::DynamicVector<double> UKinetic(const blaze::DynamicVector<double>& Lambda);
    blaze::DynamicMatrix<double> UKinetic(const blaze::DynamicMatrix<double>& Lambda);
    /**
     * calculate derivative of solution for kinetic entropy with given entropic variable
     * @param entropic variable
     * @return derivative of solution
     */
    blaze::DynamicMatrix<double> DUKinetic(const blaze::DynamicVector<double> &Lambda);
    const std::vector<blaze::DynamicVector<double>>& GetPhi(){return _phi;}
    const blaze::DynamicVector<double>& GetPhiTilde(int k){return _phiTilde[k];}
    const std::vector<blaze::DynamicVector<double>>& GetPhiTilde(){return _phiTilde;}
    /**
     * Transform matrix to vector
     * @return
     */
    blaze::DynamicVector<double> MakeVector(const blaze::DynamicMatrix<double> &mat)const;
    /**
     * Transform vector to matrix
     * @return
     */
    blaze::DynamicMatrix<double> MakeMatrix(const blaze::DynamicVector<double> &vec)const;

    double GetUPlus()const {return _uPlus;}
    double GetUMinus()const {return _uMinus;}
    double CalcNorm(blaze::DynamicVector<double> &test);
};

#endif // CLOSURE_H
