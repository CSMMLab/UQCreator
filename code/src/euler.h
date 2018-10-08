#ifndef EULER_H
#define EULER_H

#include "problem.h"

class Euler : public Problem
{
  private:
    double _gamma;

  public:
    Euler( std::string inputFile );
    virtual void Solve();
    virtual void Plot() const;
    virtual void Print();
    virtual void WriteToFile( std::string filename, int filetype ) const;
    virtual Vector G( const Vector& u, const Vector& v );
    virtual Matrix G( const Matrix& u, const Matrix& v );
    virtual void Plot( Vector& x, Vector& u );
    virtual double ExactSolution( double t, double x, double xi );
    Vector F( const Vector& u );

    Matrix F( const Matrix& u );
    double GetGamma() const;
};

#endif    // EULER_H
