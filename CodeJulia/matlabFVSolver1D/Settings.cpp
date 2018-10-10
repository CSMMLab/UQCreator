#include "Settings.h"

Settings::Settings() : _testCase( 0 ) {}

Settings::~Settings() {}

void Settings::SetDt( double dt ) { _dt = dt; }
double Settings::GetDt() { return _dt; }

void Settings::SetNTimesteps( unsigned int nTimesteps ) { _nTimesteps = nTimesteps; }

unsigned int Settings::GetNTimesteps() { return _nTimesteps; }

void Settings::SetMaxCflNumber( double maxCflNumber ) { _maxCflNumber = maxCflNumber; }

double Settings::GetMaxCflNumber() { return _maxCflNumber; }

void Settings::SetTimeIntegrator( unsigned int timeIntegratorFlag ) { _timeIntegratorFlag = timeIntegratorFlag; }

unsigned int Settings::GetTimeIntegrator() { return _timeIntegratorFlag; }

void Settings::SetLimiter( unsigned int limiterFlag ) { _limiterFlag = limiterFlag; }

unsigned int Settings::GetLimiter() { return _limiterFlag; }

void Settings::SetNOutputs( unsigned int nOutputs ) { _nOutputs = nOutputs; }

unsigned int Settings::GetNOutputs() { return _nOutputs; }

void Settings::SetOutputStepSize( unsigned int outputStepSize ) { _outputStepSize = outputStepSize; }

unsigned int Settings::GetOutputStepSize() { return _outputStepSize; }

void Settings::SetScreenOutputStepSize( unsigned int screenOutputStepSize ) { _screenOutputStepSize = screenOutputStepSize; }

unsigned int Settings::GetScreenOutputStepSize() { return _screenOutputStepSize; }

void Settings::SetGridDimension( unsigned int gridDimension ) { _gridDimension = gridDimension; }

unsigned int Settings::GetGridDimension() { return _gridDimension; }

void Settings::SetPhysicalProblem( unsigned int physicalProblemFlag ) { _physicalProblemFlag = physicalProblemFlag; }

unsigned int Settings::GetPhysicalProblem() { return _physicalProblemFlag; }

void Settings::SetNCells( unsigned int nCells ) { _nCells = nCells; }

unsigned int Settings::GetNCells() { return _nCells; }

void Settings::SetDimX( unsigned int dimX ) { _dimX = dimX; }

unsigned int Settings::GetDimX() { return _dimX; }

void Settings::SetXLeft( float xLeft ) { _xLeft = xLeft; }

float Settings::GetXLeft() { return _xLeft; }

void Settings::SetXRight( float xRight ) { _xRight = xRight; }

float Settings::GetXRight() { return _xRight; }

void Settings::SetDimY( unsigned int dimY ) { _dimY = dimY; }

unsigned int Settings::GetDimY() { return _dimY; }

void Settings::SetYLeft( float yLeft ) { _yLeft = yLeft; }

float Settings::GetYLeft() { return _yLeft; }

void Settings::SetYRight( float yRight ) { _yRight = yRight; }

float Settings::GetYRight() { return _yRight; }

void Settings::HandleHeating( bool set ) { _useHeating = set; }

bool Settings::HeatingOn() { return _useHeating; }

void Settings::SetQ0( double Q ) { _Q = Q; }

double Settings::GetQ0() { return _Q; }

void Settings::SetMaxHeatingWidth( double maxWidth ) { _maxWidth = maxWidth; }

double Settings::GetMaxHeatingWidth() { return _maxWidth; }

void Settings::SetTestCase( int tc ) { _testCase = tc; }

int Settings::GetTestCase() const { return _testCase; }

void Settings::SetBounCon( int bc ) { _bounCon = bc; }

int Settings::GetBounCon() const { return _bounCon; }

void Settings::SetScaling( bool useScaling ) { _useScaling = useScaling; }

bool Settings::ScalingOn() { return _useScaling; }
