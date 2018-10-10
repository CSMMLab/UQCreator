/** Class defines physical problem as Burger's equation
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 06.05.2013
 *
 *
 *  @version 2
 *  CFL number added
 *  @version 1
 *  Burger's equation(first definition)
 */
#ifndef BURGERSEQUATION_H
#define BURGERSEQUATION_H

#include "PhysicalProblem.h"

class BurgersEquation : public PhysicalProblem
{
  private:
    const unsigned int _StateDim = 1;    ///< dimension of problem (number of conservation equations)
  public:
    BurgersEquation( Settings* settings );    ///< constructor
    virtual ~BurgersEquation();               ///< destructor
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
    void t1_EstimateMinSpeed( Vector* state, Vector* t1_state, double& minSpeed, double& t1_minSpeed ) {}
    /**  estimator for max speed
     *  @param[in]  state vector
     *  @return value for max speed estimation
     */
    void EstimateMaxSpeed( Vector* state, double& maxSpeed );
    void t1_EstimateMaxSpeed( Vector* state, Vector* t1_state, double& maxSpeed, double& t1_maxSpeed ) {}
    /**  getter for problem dimension
     *  @return dimension of state vector
     */
    unsigned int GetStateDim() const;
    double CFL( Vector* state, double dt, double dx );    ///< calculate CFL number
};

#endif    // BURGERSEQUATION_H
