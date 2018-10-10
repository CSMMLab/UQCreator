/** class contains settings for output, numerical methods, physics, mesh
 *
 *  @author Jonas Kusch
 *  @author Jannick Wolters
 *  @author Fabian Key
 *  @date 15.05.2013
 *
 *  @version 2
 *  added defines for setting flags
 *  @version 1
 *  general settings for output, numerical methods, physics, mesh(first definition)
 */

#ifndef SETTINGS_H
#define SETTINGS_H

// defines

// limiter
#define NOLIMITER 0
#define VANLEER 1
#define MINMOD 2

// physical problem
#define EULER1DREALFLUID 0
#define EULER1DIDGAS 1
#define BURGERS 2
#define ADVECTION 3

#define OPEN 0
#define PERIODIC 1

// timeintegrator
#define HEUN 0
#define RADAU 1

class Settings
{
  private:
    // numerical settings
    unsigned int _nTimesteps;
    double _dt;
    double _maxCflNumber;
    unsigned int _timeIntegratorFlag;
    unsigned int _limiterFlag;
    // output settings
    unsigned int _nOutputs;
    unsigned int _outputStepSize;
    unsigned int _screenOutputStepSize;
    bool _useScaling;

    // physical settings
    unsigned int _gridDimension;
    unsigned int _physicalProblemFlag;
    bool _useHeating;
    double _Q;
    double _maxWidth;
    int _testCase;
    int _bounCon;

    // mesh settings
    unsigned int _nCells;
    unsigned int _dimX;
    float _xLeft, _xRight;
    unsigned int _dimY;
    float _yLeft, _yRight;

  public:
    Settings();             ///< constructor
    virtual ~Settings();    ///< destructor

    // numerical settings
    void SetNTimesteps( unsigned int nTimesteps );                ///< setter for number of timesteps
    unsigned int GetNTimesteps();                                 ///< getter for number timesteps
    void SetDt( double dt );                                      ///< setter for timestep size
    double GetDt();                                               ///< getter for timestep size
    void SetMaxCflNumber( double maxCflNumber );                  ///< setter for max CFL number
    double GetMaxCflNumber();                                     ///< getter for max CFL number
    void SetTimeIntegrator( unsigned int timeIntegratorFlag );    ///< setter for time integrator
    unsigned int GetTimeIntegrator();                             ///< getter for time integrator
    void SetLimiter( unsigned int limiterFlag );                  ///< setter for limiter
    unsigned int GetLimiter();                                    ///< getter for limiter

    // output settings
    void SetNOutputs( unsigned int nOutputs );                      ///< setter for number of outputs
    unsigned int GetNOutputs();                                     ///< getter for number of outputs
    void SetOutputStepSize( unsigned int outputStepSize );          ///< setter for outputStepSize
    unsigned int GetOutputStepSize();                               ///< getter for outputStepSize
    void SetScreenOutputStepSize( unsigned int outputStepSize );    ///< setter for outputStepSize
    unsigned int GetScreenOutputStepSize();                         ///< getter for outputStepSize
    void SetScaling( bool useScaling );                             ///< setter for adaptive plotsize
    bool ScalingOn();                                               ///< getter fpr adaptive plotsize

    // physical settings
    void SetGridDimension( unsigned int gridDimension );            ///< setter for grid dimension
    unsigned int GetGridDimension();                                ///< getter for grid dimension
    void SetPhysicalProblem( unsigned int physicalProblemFlag );    ///< setter for physical problem flag
    unsigned int GetPhysicalProblem();                              ///< getter for physical problem flag
    void HandleHeating( bool set );                                 ///< setter for physical problem flag
    bool HeatingOn();                                               ///< getter for physical problem flag
    void SetQ0( double Q );
    double GetQ0();
    void SetMaxHeatingWidth( double maxWidth );
    double GetMaxHeatingWidth();
    void SetTestCase( int tc );
    int GetTestCase() const;
    void SetBounCon( int bc );
    int GetBounCon() const;

    // mesh settings
    void SetNCells( unsigned int nCells );    ///< setter for number of cells
    unsigned int GetNCells();                 ///< getter for number of cells
    void SetDimX( unsigned int nCells );      ///< setter for number of cells in x-direction
    unsigned int GetDimX();                   ///< getter for number of cells in x-direction
    void SetXLeft( float nCells );            ///< setter for left border in x-direction
    float GetXLeft();                         ///< getter for left border in x-direction
    void SetXRight( float nCells );           ///< setter for right border in x-direction
    float GetXRight();                        ///< getter for right border in x-direction
    void SetDimY( unsigned int nCells );      ///< setter for number of cells in y-direction
    unsigned int GetDimY();                   ///< getter for number of cells in y-direction
    void SetYLeft( float nCells );            ///< setter for left border in y-direction
    float GetYLeft();                         ///< getter for left border in y-direction
    void SetYRight( float nCells );           ///< setter for right border in y-direction
    float GetYRight();                        ///< getter for right border in y-direction
};

#endif    // SETTINGS_H
