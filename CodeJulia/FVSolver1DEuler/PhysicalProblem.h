/** Class defines physical problem
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 25.04.2013
 *
 *  @version 2
 *  added pointer of class 'Settings'
 *  @version 1
 *  definition of physical problem (first definition)
 */

#ifndef PHYSICALPROBLEM_H
#define PHYSICALPROBLEM_H

#include "Settings.h"
#include "Vector.h"

class PhysicalProblem
{
  protected:
    Settings* _settings;    ///< grants access to settings
  private:
    const unsigned int _StateDim = 1;    ///< dimension of problem (number of conservation equations)
    /** primitive variables
     * @li pressure
     * @li temperature
     * @li velocity in x-direction
     * @li density
     * @li speed of sound
     */
    double _p, _u, _rho, _a;
    /**  transformation from conservative to primitive variables
     *  @param[in]  state vector
     */
    virtual void TransformToPrimitives( Vector* state );
    virtual void t1_TransformToPrimitives( Vector* state, Vector* t1_state );

  public:
    PhysicalProblem( Settings* settings );    ///< constructor
    virtual ~PhysicalProblem();               ///< destructor
    /**  analytic flux is calculated
     *  @param[in]  state vector
     *  @param[out] vector for analytic flux evaluated at state
     */
    virtual void AnalyticFlux( Vector* state, Vector* flux )                                       = 0;
    virtual void t1_AnalyticFlux( Vector* state, Vector* t1_state, Vector* flux, Vector* t1_flux ) = 0;
    /**  estimator for min speed
     *  @param[in]  state vector
     *  @return value for min speed estimation
     */
    virtual void EstimateMinSpeed( Vector* state, double& minSpeed )                                           = 0;
    virtual void t1_EstimateMinSpeed( Vector* state, Vector* t1_state, double& minSpeed, double& t1_minSpeed ) = 0;
    /**  estimator for max speed
     *  @param[in]  state vector
     *  @return value for max speed estimation
     */
    virtual void EstimateMaxSpeed( Vector* state, double& maxSpeed )                                           = 0;
    virtual void t1_EstimateMaxSpeed( Vector* state, Vector* t1_state, double& maxSpeed, double& t1_maxSpeed ) = 0;
    /**  calculate CFL number
     *  @param[in]  state vector
     *  @param[in]  timestep size
     *  @param[in]  cell width
     */
    virtual double CFL( Vector* state, double dt, double dx ) = 0;
    /**  getter for problem dimension
     *  @return dimension of state vector
     */
    virtual unsigned int GetStateDim() const = 0;
    /**  transformation from primitive to conservative variables
     *  @param[in]  density
     *  @param[in]  pressure
     *  @param[in]  velocity
     *  @param[out]  state vector
     */
    virtual void TransformToConservatives( double rho, double p, double u, Vector* state );
    /**
     *  calculation of primitive variables
     *  @param[in]  state vector
     *  @param[out]  pressure
     *  @param[out]  velocity
     *
     */
    virtual void GetPrimitives( Vector* state, double& p, double& u );
    /**
     *  calculation of different variables
     *  @param[in]  state vector
     *  @param[out]  Temperature
     *  @param[out]  Mach number
     *  @param[out]  vapour fraction
     *
     */
    virtual void GetMisc( Vector* state, double& T, double& mach, double& vapour );
    /**
     *  calculation of different variables
     *  @param[in]  state vector
     *  @param[out]  Temperature
     *  @param[out]  Mach number
     *
     */
    virtual void GetMisc( Vector* state, double& T, double& mach );
    virtual double GetPressure( Vector* state );          ///< calculate pressure
    virtual double GetVelocity( Vector* state );          ///< calculate velocity
    virtual double GetTemperature( Vector* state );       ///< calculate temperature
    virtual double GetMachNumber( Vector* state );        ///< calculate Mach number
    virtual double GetVapourFraction( Vector* state );    ///< calculate vapour fraction
    /**
     * @brief Function gives HeatFlux, which is applied to the tube.
     *        The heat flux is simply added to the right hand side. It is not jet implemented for Advection and Burgers.
     *        It is only used if Rhs1DHeating is used.
     * @param[in] x Coordinate of cell
     * @return value for heat flux
     */
    virtual double HeatFlux( double x ) const;
    double _cfl;
};

#endif    // PHYSICALPROBLEM_H
