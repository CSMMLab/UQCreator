/** Abstract class for Rieamann solver (returns numerical flux over one interface)
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *
 *  @version 1
 *  abstract Riemann solver (first definition)
 */

#ifndef RIEMANNSOLVER_H
#define RIEMANNSOLVER_H

#include "Vector.h"

class RiemannSolver
{
  private:
  public:
    RiemannSolver();             ///< constructor
    virtual ~RiemannSolver();    ///< destructor

    /** calculates the numerical flux between two cells
     * @param[in] left, state vector of the cell left to the Riemann problem
     * @param[in] right, state vector of the cell right to the Riemann problem
     * @param[out] numericalFlux, solution vector (hand in by reference)
     */
    virtual void NumericalFlux( Vector* left, Vector* right, Vector* numericalFlux ) = 0;
    virtual void
    t1_NumericalFlux( Vector* left, Vector* t1_left, Vector* right, Vector* t1_right, Vector* numericalFlux, Vector* t1_numericalFlux ) = 0;
};

#endif    // RIEMANNSOLVER_H
