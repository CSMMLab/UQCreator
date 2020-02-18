#ifndef KINETICCLOSURE_H
#define KINETICCLOSURE_H

#include "closure.h"

class KineticClosure : public Closure
{
  private:
    KineticClosure() = delete;
    unsigned _nMoments;

  public:
    KineticClosure( Settings* settings );
    virtual ~KineticClosure();

    virtual void U( Vector& out, const Vector& Lambda );
    virtual void U( Vector& out, const Vector& Lambda, bool dummy );
    virtual void U( Matrix& out, const Matrix& Lambda );
    virtual Matrix U( const Matrix& Lambda );
    virtual void DU( Matrix& y, const Vector& Lambda );
};

#endif // KINETICCLOSURE_H
