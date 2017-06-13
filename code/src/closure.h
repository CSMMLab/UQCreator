#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <vector>

typedef blaze::DynamicVector<double> vector;
typedef blaze::DynamicMatrix<double> matrix;

class Closure
{
private:
    Problem* _problem;
    Quadrature* _quadrature;
    Newton* _newton;
    BasisFunctions* _basis;
    std::vector<vector> _phi; // stores basis functions evaluated at quadrature points
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
     * @return The test results
     */
    Closure(Problem *problem);
    /**
     * calculate dual vector fulfilling the moment constraint
     * @param moment vector
     * @param initial guess for dual vector
     * @return correct dual vector
     */
    vector SolveClosure(vector u, vector lambda);
    double EvaluateLambda(vector lambda,double xi);
    double UKinetic(double Lambda);
    double DUKinetic(double Lambda);
};

#endif // CLOSURE_H
