/** Abstract class for limiter
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 06.05.2013
 *
 *
 *  @version 1
 *  abstract limiter(first definition)
 */

#ifndef LIMITER_H
#define LIMITER_H

#include "Mesh.h"
#include "PhysicalProblem.h"
#include "Vector.h"

class Limiter
{
  private:
  protected:
    PhysicalProblem* _physicalProblem;

  public:
    Limiter( PhysicalProblem* physicalProblem );    ///< constructor
    virtual ~Limiter();                             ///< destructor
    /**
     * calculates limited state values at left and right interface
     * @param[in] index, cell index
     * @param[in] mesh, mesh
     * @param[out] limitLeft, limited state at left interface
     * @param[out] limitRight, limited state at right interface
     */
    virtual void CalculateLimitedState1D( int index, Mesh* mesh, Vector* limitLeft, Vector* limitRight ) = 0;
    virtual void t1_CalculateLimitedState1D(
        int index, Mesh* mesh, Mesh* t1_mesh, Vector* limitLeft, Vector* t1_limitLeft, Vector* limitRight, Vector* t1_limitRight ) = 0;
};
#endif    // LIMITER_H
