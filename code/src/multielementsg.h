#ifndef MULTIELEMENTSG_H
#define MULTIELEMENTSG_H

#include "closure.h"

class MultiElementSG : public Closure
{
private:
  MultiElementSG() = delete;
public:
  MultiElementSG( Settings* settings );
  virtual ~MultiElementSG();

  virtual void U( Vector& out, const Vector& Lambda );
  virtual void U( Matrix& out, const Matrix& Lambda );
  virtual Matrix U( const Matrix& Lambda );
  virtual void DU( Matrix& y, const Vector& Lambda );
  virtual void SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel );
  virtual void SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel );
  Vector EvaluateLambda( const Matrix& lambda, unsigned k, unsigned nTotal );
  Matrix EvaluateLambda( const Matrix& lambda ) const;
  void EvaluateLambda( Matrix& out, const Matrix& lambda ) const;
  Matrix EvaluateLambdaOnPE( const Matrix& lambda, unsigned levelOld, unsigned levelNew ) const;
};

#endif // MULTIELEMENTSG_H
