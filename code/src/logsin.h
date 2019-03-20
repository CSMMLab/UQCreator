#ifndef LOGSIN_H
#define LOGSIN_H

#include "closure.h"
#include "problem.h"

class LogSin : public Closure
{
  private:
    double _uMinus, _uPlus;
    LogSin() = delete;

  public:
    LogSin( Settings* settings );
    virtual ~LogSin();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // LOGSIN_H
