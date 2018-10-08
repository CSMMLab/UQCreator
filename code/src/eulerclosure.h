#ifndef EULERCLOSURE_H
#define EULERCLOSURE_H

#include "closure.h"
#include "problem.h"

class EulerClosure : public Closure
{
private:
    double _gamma;

public:
    EulerClosure() = delete;
    EulerClosure(Problem* problem);
    virtual ~EulerClosure();

    virtual void U( blaze::DynamicVector<double>& out, const blaze::DynamicVector<double>& Lambda );
    virtual blaze::DynamicMatrix<double> U( const blaze::DynamicMatrix<double>& Lambda );
    virtual void DU( blaze::DynamicMatrix<double>& y, const blaze::DynamicVector<double>& Lambda );
};

#endif // EULERCLOSURE_H
