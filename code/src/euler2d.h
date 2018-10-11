#ifndef EULER2D_H
#define EULER2D_H

#include "problem.h"

class Euler2D : public Problem
{
  private:
    double _gamma;

  public:
    Euler2D( std::string inputFile );
    virtual void Solve();
    virtual void Print();
    virtual void WriteToFile( std::string filename, int filetype ) const;
    virtual Vector G(const Vector& u, const Vector& v , const Vector &nUnit, const Vector &n);
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    virtual double ExactSolution( double t, double x, double xi );
    Matrix F( const Vector& u );

    Matrix F( const Matrix& u );
    double GetGamma() const;
};

#endif // EULER2D_H
