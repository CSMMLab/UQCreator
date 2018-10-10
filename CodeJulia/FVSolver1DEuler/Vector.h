/** Class realizes a vector
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 28.04.2013
 *
 *
 *  @version 1
 *  vector class for saving state (first definition)
 */

#ifndef vector_inc
#define vector_inc

#include <assert.h>
#include <cmath>
#include <iostream>

class Vector
{
  private:
    /**
     * stores values
     */
    double* _a;
    /**
     * stores dimensions
     */
    unsigned _n;

  public:
    /**
     * constructor
     * @param[in] N, N is dimension of Vector
     */
    Vector( const unsigned N );
    /**
     * destructor
     */
    ~Vector();
    /**
     * copy constructor
     * @param[in] R, R is Vector which the copy constructor copies
     */
    Vector( const Vector& R );
    /**
     * = operator
     * @param[in] R, R is vector for = operator
     * @return
     */
    Vector& operator=( const Vector& R );
    /**
     * get and set [] operator
     * @param[in] i, i is index in vector
     * @return returns value at position i
     */
    double& operator[]( const unsigned i );
    const double& operator[]( const unsigned i ) const;
    /**
     * Vector-Vector operators +
     * @param[in] B, B is Vector for operation of Vectoraddition, -substitution
     * @return resulting vector from Vector operation
     */
    Vector operator+( const Vector& B ) const;
    /**
     * Vector-Vector operators -
     * @param[in] B, B is Vector for operation of Vectoraddition, -substitution
     * @return resulting vector from Vector operation
     */
    Vector operator-( const Vector& B ) const;
    /**
     * Vector-Scalar operator *
     * @param[in] Scalar, scalar for vector multiplication
     * @return resulting vector from Vector-scalar operation
     */
    Vector operator*( const double Scalar ) const;
    /**
     * returns the vector dimension
     * @return vector dimension
     */
    unsigned Dimension() const { return _n; }
    /**
     * prints the vector on console
     */
    void Show() const;
    /**
     * Addition of this vector with vector a in a more efficient way
     * @param[in] a, second vector for addition
     */
    void Add( Vector* a );
    /**
     * Subtraction of this vector with vector a in a more efficient way
     * @param[in] a, second vector for subtraction
     */
    void Subtract( Vector* a );
    /**
     * Vector-Scalar multiplication of this vector with input scalar
     * @param[in] Scalar, scalar value for multiplication
     */
    void Multiply( const double Scalar );
    /**
     * sets at values in the vector to zero
     */
    void MakeZero();

  private:
    Vector() {}
};

#endif
