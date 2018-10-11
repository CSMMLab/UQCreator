#ifndef EULER2D_H
#define EULER2D_H

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
    virtual Matrix G( const Matrix& u, const Matrix& v );
    virtual double ExactSolution( double t, double x, double xi );
    Vector F( const Vector& u );

    Matrix F( const Matrix& u );
    double GetGamma() const;
};

#endif // EULER2D_H
