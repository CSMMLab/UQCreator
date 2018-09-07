#ifndef BURGERS_H
#define BURGERS_H

#include <fstream>
#include <iostream>

#include "gnuplot-iostream.h"
#include "problem.h"

class Burgers : public Problem
{
  private:
    blaze::DynamicVector<double> _u;
    blaze::DynamicVector<double> _x;
    double _dx;
    double _dt;
    int _nCells;
    int _nTimeSteps;
    double H( double u, double v, double w );
    double F( double u );
    blaze::DynamicMatrix<double> F( const blaze::DynamicMatrix<double>& u );
    double IC( double x, double uL, double uR );

    Burgers() {}

  public:
    Burgers( std::string inputFile );
    virtual void Solve();
    virtual void Plot( blaze::DynamicVector<double>& x, blaze::DynamicVector<double>& u );
    virtual void Print();
    virtual void WriteToFile( std::string filename, int filetype ) const;
    double G( double u, double v );
    blaze::DynamicMatrix<double> G( const blaze::DynamicMatrix<double>& u, const blaze::DynamicMatrix<double>& v );
    virtual double ExactSolution( double t, double x, double xi );
};

#endif    // BURGERS_H
