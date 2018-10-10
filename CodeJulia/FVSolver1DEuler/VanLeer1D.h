/** van Leer limiter for 1D case
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 06.05.2013
 *
 *
 *  @version 1
 *  van Leer limiter(first definition)
 */

#ifndef VANLEER1D_H
#define VANLEER1D_H

#include "Limiter.h"

class VanLeer1D : public Limiter
{
  private:
  public:
    VanLeer1D( PhysicalProblem* physicalProblem );    ///< constructor
    virtual ~VanLeer1D();                             ///< destructor
    /**
     * calculates limited state values at left and right interface of one cell
     * @param[in] index, cell index
     * @param[in] mesh, mesh
     * @param[out] limitLeft, limited state at left interface
     * @param[out] limitRight, limited state at right interface
     */
    virtual void CalculateLimitedState1D( int index, Mesh* mesh, Vector* limitLeft, Vector* limitRight );
    virtual void t1_CalculateLimitedState1D(
        int index, Mesh* mesh, Mesh* t1_mesh, Vector* limitLeft, Vector* t1_limitLeft, Vector* limitRight, Vector* t1_limitRight );
};
#endif    // VANLEER1D_H
