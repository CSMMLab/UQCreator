/** Class implements HLL Riemann solver
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *
 *  @version 2
 *  function 'NumericalFlux' is now faster and more efficient, results have even been verified by hand
 *
 *  @version 1
 *  HLL (Harten, Lax, van Leer) Riemann solver (first definition)
 */
#ifndef HLLSOLVER1D_H
#define HLLSOLVER1D_H

#include "PhysicalProblem.h"
#include "RiemannSolver.h"

class HLLSolver1D : public RiemannSolver
{
  private:
    PhysicalProblem* _euler1D;    ///< allows easy access to the physical problem (handed in via construktor)
  public:
    /** constructor
     * @param[in] euler1D, physical problem for which the Riemann problem shall be solved
     */
    HLLSolver1D( PhysicalProblem* euler1D );

    virtual ~HLLSolver1D();    ///< destructor

    /** calculates the numerical flux between two cells
     * @param[in] left, state vector of the cell left to the Riemann problem
     * @param[in] right, state vector of the cell right to the Riemann problem
     * @param[out] numericalFlux, solution vector (hand in by reference)
     */
    void NumericalFlux( Vector* left, Vector* right, Vector* numericalFlux );
    void t1_NumericalFlux( Vector* left, Vector* t1_left, Vector* right, Vector* t1_right, Vector* numericalFlux, Vector* t1_numericalFlux );

  private:
    HLLSolver1D() {}    ///< invalid constructor
};

#endif    // HLLSOLVER1D_H
