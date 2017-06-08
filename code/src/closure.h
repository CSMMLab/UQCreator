#ifndef CLOSURE_H
#define CLOSURE_H

#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#include <vector>


class Closure
{
private:
    Problem* _problem;
    Quadrature* _quadrature;
    Newton* _newton;
    BasisFunctions* _basis;
    blaze::DynamicMatrix<float,blaze::rowMajor>* _phi; // stores basis functions evaluated at quadrature points
    std::vector<blaze::DynamicMatrix<float,blaze::rowMajor>*> _hPartial; // stores partial matrices for Hessian computation
public:
    Closure(Problem problem);
    blaze::DynamicVector<double,blaze::rowVector> SolveClosure(blaze::DynamicVector<double,blaze::rowVector> u);
    double uKinetic(double Lambda);

};

#endif // CLOSURE_H
