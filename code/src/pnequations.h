#ifndef PNEQUATIONS_H
#define PNEQUATIONS_H

#include "problem.h"

class PNEquations : public Problem
{
  private:
    // moment orders for P_N
    const int _N;

    //System Matrix for x and y flux
    Matrix _Ax;
    Matrix _Ay;

    // parameter functions for setting up system matrix
    double AParam(int l, int k)const;
    double BParam(int l, int k)const;
    double CParam(int l, int k)const;
    double DParam(int l, int k)const;
    double EParam(int l, int k)const;
    double FParam(int l, int k)const;

    double CTilde(int l, int k)const;
    double DTilde(int l, int k)const;
    double ETilde(int l, int k)const;
    double FTilde(int l, int k)const;

    // mathematical + index functions
    int Sgn(int k)const;
    int kPlus(int k)const;
    int kMinus(int k)const;
    unsigned GlobalIndex(int l, int k)const;

    // function for setting up system matrices
    void SetupSystemMatrices();

  public:
    PNEquations( Settings* settings );
    virtual ~PNEquations();
    inline Vector G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n );
    virtual Matrix G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n );
    Matrix F( const Vector& u );
    Matrix F( const Matrix& u );
    virtual double ComputeDt( const Matrix& u, double dx ) const;
    virtual Vector IC( const Vector& x, const Vector& xi );
};

#endif // PNEQUATIONS_H
