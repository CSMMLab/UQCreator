/** Class defines physical problem as euler eqautions for ideal gases
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 06.05.2013
 *
 *
 *  @version 1
 *  euler equations for ideal gases (first definition)
 */

#ifndef EULEREQUATIONSIDGAS_H
#define EULEREQUATIONSIDGAS_H

#include "PhysicalProblem.h"

class EulerEquationsIdGas : public PhysicalProblem
{
  private:
    const unsigned int _StateDim = 3;    ///< dimension of problem (number of conservation equations)

    // const double _Gamma=1.4;///< isentropic exponent
    // const double _SpecificR=287;// J/(kg*K)
    const double _Gamma     = 1.333333333333333;    ///< isentropic exponent
    const double _SpecificR = 462;                  // J/(kg*K)
    /**
     * specifications for heat flux
     * @brief q0 maximal heating at x = 0
     * @brief maxWidth width of the maxwellian
     */
    double _q0;
    double _maxWidth;

    /** primitive variables
     * @li pressure
     * @li velocity in x-direction
     * @li density
     * @li speed of sound
     * @li temperature
     */
    double _p, _u, _rho, _a, _T;
    double _t1_p, _t1_u, _t1_rho, _t1_a;
    /**  transformation from conservative to primitive variables
     *  @param[in]  state vector
     */
    void TransformToPrimitives( Vector* state );
    void t1_TransformToPrimitives( Vector* state, Vector* t1_state ) {}

  public:
    EulerEquationsIdGas( Settings* settings );    ///< constructor
    virtual ~EulerEquationsIdGas();               ///< destructor
    /**  analytic flux is calculated
     *  @param[in]  state vector
     *  @param[out] vector for analytic flux evaluated at state
     */
    void AnalyticFlux( Vector* state, Vector* flux );
    void t1_AnalyticFlux( Vector* state, Vector* t1_state, Vector* flux, Vector* t1_flux ) {}
    /**  estimator for min speed
     *  @param[in]  state vector
     *  @return value for min speed estimation
     */
    void EstimateMinSpeed( Vector* state, double& minSpeed );
    void t1_EstimateMinSpeed( Vector* state, Vector* t1_state, double& minSpeed, double& t1_minSpeed );
    /**  estimator for max speed
     *  @param[in]  state vector
     *  @return value for max speed estimation
     */
    void EstimateMaxSpeed( Vector* state, double& maxSpeed );
    void t1_EstimateMaxSpeed( Vector* state, Vector* t1_state, double& maxSpeed, double& t1_maxSpeed );
    /**  calculate CFL number
     *  @param[in]  state vector
     *  @param[in]  timestep size
     *  @param[in]  cell width
     */
    virtual double CFL( Vector* state, double dt, double dx );
    /**  getter for problem dimension
     *  @return dimension of state vector
     */
    unsigned int GetStateDim() const;
    /**  transformation from primitive to conservative variables
     *  @param[in]  density
     *  @param[in]  pressure
     *  @param[in]  velocity
     *  @param[out]  state vector
     */
    void TransformToConservatives( double rho, double p, double u, Vector* state );
    /**
     *  calculation of primitive variables
     *  @param[in]  state vector
     *  @param[out]  pressure
     *  @param[out]  velocity
     *
     */
    void GetPrimitives( Vector* state, double& p, double& u );
    /**
     *  calculation of different variables
     *  @param[in]  state vector
     *  @param[out]  Temperature
     *  @param[out]  Mach number
     *
     */
    void GetMisc( Vector* state, double& T, double& Mach );
    double GetPressure( Vector* state );       ///< calculate pressure
    double GetVelocity( Vector* state );       ///< calculate velocity
    double GetTemperature( Vector* state );    ///< calculate temperature
    double GetMachNumber( Vector* state );     ///< calculate Mach number
    /**
     * @brief Function gives HeatFlux, which is applied to the tube.
     *        The heat flux is simply added to the right hand side. It is not jet implemented for Advection and Burgers.
     *        It is only used if Rhs1DHeating is used.
     * @param[in] x Coordinate of cell
     * @return value for heat flux
     */
    virtual double HeatFlux( double x ) const;
};

#endif    // EULEREQUATIONSIDGAS_H
