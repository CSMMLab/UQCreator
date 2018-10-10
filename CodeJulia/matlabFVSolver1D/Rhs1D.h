/** Class implements 1D right hand side for ODE system
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *  @version 3
 *  SetupRHS with index i and mesh as input added. This is easier for the limiter.
 *  @version 2
 *  comments added, SetupRhs completed, _mesh is now saved in Rhs1D as Mesh1D
 *  @version 1
 *  right hand side for ODE system 1D case (first definition)
 */

#ifndef RHS1D_H
#define RHS1D_H

#include "HLLSolver1D.h"
#include "Mesh1D.h"
#include "Rhs.h"
#include "Vector.h"

class Rhs1D : public Rhs
{
  private:
    Vector *_limitRightLeft, *_limitLeftCurrent, *_limitRightCurrent, *_limitLeftRight;
    Vector *_t1_limitRightLeft, *_t1_limitLeftCurrent, *_t1_limitRightCurrent, *_t1_limitLeftRight;

  public:
    /**
     * constructor
     * @param[in] physProb, stores physical problem
     */
    Rhs1D( Settings* settings, PhysicalProblem* physProb, double dx );
    /**
     * destructor
     */
    virtual ~Rhs1D();
    /**
     * Function for calculation of right hand side
     * @param[in] index, Calculates RHS for Cell i
     * @param[in] mesh, pointer to mesh on which cell i can be found
     * @param[out] result, pointer to Vector with values of the Right hand side
     */
    void SetupRhs( int index, Mesh* mesh, Vector* result );
    void t1_SetupRhs( int index, Mesh* mesh, Mesh* t1_mesh, Vector* result, Vector* t1_result );
};

#endif    // RHS1D_H
