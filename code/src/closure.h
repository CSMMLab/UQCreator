#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/LAPACK.h>
#include <vector>
#include "quadrature.h"
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
    blaze::DynamicMatrix<double> Hessian(blaze::DynamicVector<double> lambda);
    blaze::DynamicVector<double> Gradient(blaze::DynamicVector<double> lambda, blaze::DynamicVector<double> u);
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
    blaze::DynamicVector<double> SolveClosure(blaze::DynamicVector<double> u, blaze::DynamicVector<double> lambda);
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    double EvaluateLambda(blaze::DynamicVector<double> lambda,double xi);
    /**
     * calculate solution for kinetic entropy with given entropic variable
     * @param entropic variable
     * @return solution
     */
    double UKinetic(double Lambda);
    /**
     * calculate derivative of solution for kinetic entropy with given entropic variable
     * @param entropic variable
     * @return derivative of solution
     */
    double DUKinetic(double Lambda);
    std::vector<blaze::DynamicVector<double>> GetPhi(){return _phi;}
    blaze::DynamicVector<double> GetPhiTilde(int k){return _phiTilde[k];}
    std::vector<blaze::DynamicVector<double>> GetPhiTilde(){return _phiTilde;}
};

#endif // CLOSURE_H
