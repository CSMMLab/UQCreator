/** Class implements Heun time integrator
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *  @version 2
 *  deleted parameter 'dt', timestep size by '_settings'
 *  @version 1
 *  Heun time integrator(first definition)
 */

#ifndef HEUN_H
#define HEUN_H

#include "TimeIntegrator.h"

class Heun : public TimeIntegrator
{
  private:
    Mesh* _intermediates;
    PhysicalProblem* _physicalProblem;
    /**
     * Calculates cell intermediate x_i+1/2. Used for Integrate.
     * @param[in] i, index of cell which needs intermediate
     * @param[out] intermediate, Vector in which result is stored
     */
    void CalcIntermediate( int i, Vector* intermediate );
    /**
     * Calculates and stores time update for cell[index]
     * @param[in] index, index of cell which is to be updated
     */
    void UpdateWithIntermediate( int index );

  public:
    /**
     * Constructor for Heun Class
     * @param[in] dt, length of one time step ( Heun has constant time step )
     * @param[in] mesh, input data structure with cell values to update
     */
    Heun( Settings* settings, Mesh* mesh, PhysicalProblem* physicalProblem, Rhs* rhs );
    /**
     * destructor of class Heun
     */
    virtual ~Heun();
    /**
     * Integrates to next time step and updates cell values.
     * Uses second mesh to store intermediate values.
     */
    void Integrate();

  private:
    /**
     * invalid constructor for Heun
     */
    // Heun(): TimeIntegrator() {}
};

#endif    // HEUN_H
