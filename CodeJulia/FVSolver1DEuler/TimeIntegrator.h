/** Abstract class for time integrator
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 29.04.2013
 *
 *  @version 2
 *  added pointer of class 'Settings' which contains information about timestep size
 *  @version 1
 *  abstract time integrator (first definition)
 */

#ifndef TIMEINTEGRATOR_H
#define TIMEINTEGRATOR_H

#include "Mesh.h"
#include "Rhs.h"
#include "Settings.h"

class TimeIntegrator
{
  private:
  protected:
    Settings* _settings;    ///< grants access to settings
    Mesh* _mesh;            ///< grants access to the mesh
    Rhs* _rhs;              ///< grants access to the right hand side
  public:
    TimeIntegrator( Settings* settings, Mesh* mesh, Rhs* rhs );    ///< constructor
    virtual ~TimeIntegrator();                                     ///< destructor
    virtual void Integrate() = 0;    /// abstract function which is supposed to execute one time step along the whole mesh
};

#endif    // TIMEINTEGRATOR_H
