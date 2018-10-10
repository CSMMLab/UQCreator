/** Class implements abstract right hand side for ODE system
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *
 *  @version 1
 *  abstract right hand side for ODE system (first definition)
 */

#ifndef RHS_H
#define RHS_H

#include "Limiter.h"
#include "Mesh.h"
#include "PhysicalProblem.h"
#include "RiemannSolver.h"

class Rhs
{
  protected:
    Settings* _settings;
    /**
     * stores Riemann Solver for calculation of Right Hand Side
     */
    RiemannSolver* _riemannSolver;
    /**
     * stores Limiter for calculation of Right Hand Side
     */
    Limiter* _limiter;
    /**
     * stores Physical Problem especially for calling the Riemann solver
     */
    PhysicalProblem* _physProb;
    /**
     * stores step-/cellsize
     */
    double _dx;

  public:
    /**
     * constructor
     * @param[in] mesh, geometry containing cell values which are needed for right hand side
     * @param[in] physProb, stores physical problem
     */
    Rhs( Settings* settings, PhysicalProblem* physProb, double dx );
    /**
     * destructor for Rhs
     */
    virtual ~Rhs();
    /**
     * Abstract function for calculation of right hand side
     * @param[in] index, Calculates RHS for Cell i
     * @param[in] mesh, pointer to mesh on which cell i can be found
     * @param[out] result, pointer to Vector with values of the Right hand side
     */
    virtual void SetupRhs( int index, Mesh* mesh, Vector* result )                                      = 0;
    virtual void t1_SetupRhs( int index, Mesh* mesh, Mesh* t1_mesh, Vector* result, Vector* t1_result ) = 0;

  private:
    /**
     * Invalid constructor
     */
    Rhs() {}
};

#endif    // RHS_H
