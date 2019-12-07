#ifndef LOGBARRIERCLOSURE_H
#define LOGBARRIERCLOSURE_H

#include "closure.h"
#include "problem.h"

class LogBarrierClosure : public Closure
{
  private:
    double _uMinus, _uPlus;
    LogBarrierClosure() = delete;

    double U( const double Lambda );

  public:
    LogBarrierClosure( Settings* settings );
    virtual ~LogBarrierClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );

    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // LOGBARRIERCLOSURE_H


