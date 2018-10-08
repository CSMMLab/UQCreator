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

    virtual void U( Vector& out, const Vector& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // EULERCLOSURE_H
