#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <vector>
#include "quadrature.h"
#include "legendre.h"

typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;

class Closure
{
private:
    Problem* _problem;
    Quadrature* _quadrature;
    Polynomial* _basis;
    std::vector<vector> _phi; // stores basis functions evaluated at quadrature points
    std::vector<vector> _phiTilde; // stores scaled basis functions evaluated at quadrature points
    std::vector<matrix> _hPartial; // stores partial matrices for Hessian computation
    double _uMinus, _uPlus; // IPM bounds for scalar problems
    int _nMoments;
    int _nQuadPoints;
    matrix Hessian(vector lambda);
    vector Gradient(vector lambda, vector u);
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
    vector SolveClosure(vector u, vector lambda);
    /**
     * calculate entropic variable from given dual vector
     * @param dual variable
     * @param dual variable
     * @return entropic state
     */
    double EvaluateLambda(vector lambda,double xi);
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
    std::vector<vector> GetPhi(){return _phi;}
    vector GetPhiTilde(int k){return _phiTilde[k];}
};

#endif // CLOSURE_H
