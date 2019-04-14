#ifndef REGULARIZEDEULER_H
#define REGULARIZEDEULER_H

#include "closure.h"

class RegularizedEuler : public Closure
{
private:
  double _gamma;
  double _eta; // regularization parameter
  RegularizedEuler() = delete;
  virtual void Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal );
  virtual void Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal );
public:
  RegularizedEuler( Settings* settings );
  virtual ~RegularizedEuler();

  virtual void U( Vector& out, const Vector& Lambda );
  virtual void U( Matrix& out, const Matrix& Lambda );
  virtual Matrix U( const Matrix& Lambda );
  virtual void DU( Matrix& y, const Vector& Lambda );
  virtual void DS( Vector& ds, const Vector& u ) const;
};

#endif // REGULARIZEDEULER_H
