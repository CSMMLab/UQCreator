#define N_prob 0
#define N_init 1
#define N_dt 2
#define N_nSteps 3
#define N_maxCFl 4
#define N_BC 5
#define N_Limiter 6
#define N_bounds 7

#include "Advection.h"
#include "BurgersEquation.h"
#include "EulerEquationsIdGas.h"
#include "Heun.h"
#include "Mesh1D.h"
#include "Rhs1D.h"
#include "matrix.h"
#include "mex.h"

#include <cstring>
#include <iostream>
#include <string>

// these casts are dirty af
double getScalar( int n, const mxArray* prhs[] ) { return ( (double)*mxGetPr( prhs[n] ) ); }

double* getArray( int n, int dim, const mxArray* prhs[] ) {
    auto tmp = const_cast<mxArray**>( prhs );
    mxArray* lhs[1];
    mexCallMATLAB( 1, lhs, 1, &( tmp[n] ), "double" );
    double* ret = (double*)mxGetData( lhs[0] );
    return ret;
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] ) {
    int nStates = mxGetN( prhs[N_init] );
    int nCells  = mxGetM( prhs[N_init] );

    int n1 = mxGetN( prhs[N_bounds] );
    int n2 = mxGetM( prhs[N_bounds] );

    assert( nStates > 0 );
    double* ics = getArray( N_init, nCells * nStates, prhs );
    std::vector<std::vector<double>> _initialCellStates( nCells, std::vector<double>( nStates, 0.0 ) );
    for( int i = 0; i < nStates; ++i ) {
        for( int j = 0; j < nCells; ++j ) {
            _initialCellStates[j].at( i ) = ics[j + i * nCells];
        }
    }

    double maxCFL, dt, nSteps, BC, Limiter, problem;
    double* bounds;
    problem = getScalar( N_prob, prhs );
    dt      = getScalar( N_dt, prhs );
    nSteps  = getScalar( N_nSteps, prhs );
    maxCFL  = getScalar( N_maxCFl, prhs );
    BC      = getScalar( N_BC, prhs );
    Limiter = getScalar( N_Limiter, prhs );
    bounds  = getArray( N_bounds, 2, prhs );
    /*
        mexPrintf( "%f\n", problem );
        mexPrintf( "%f\n", dt );
        mexPrintf( "%f\n", nSteps );
        mexPrintf( "%f\n", maxCFL );
        mexPrintf( "%f\n", BC );
        mexPrintf( "%f\n", Limiter );
        mexPrintf( "%f\n", bounds[0] );
        mexPrintf( "%f\n", bounds[1] );
    */
    ///////////
    // inits //
    ///////////

    Settings* _settings = new Settings;
    TimeIntegrator* _timeIntegrator;
    Mesh* _mesh;
    Rhs* _rhs;
    PhysicalProblem* _physicalProblem;

    // hardcoded settings for matlab use
    _settings->SetBounCon( BC );
    _settings->SetGridDimension( 1 );
    _settings->SetDimX( nCells );
    _settings->SetXLeft( bounds[0] );
    _settings->SetXRight( bounds[1] );
    _settings->SetDimY( 0 );
    _settings->SetYLeft( 0 );
    _settings->SetYRight( 0 );
    _settings->SetDt( dt );
    _settings->SetLimiter( Limiter );
    _settings->SetMaxCflNumber( maxCFL );
    _settings->SetNTimesteps( nSteps );
    _settings->HandleHeating( false );
    _settings->SetQ0( 0.0 );
    _settings->SetMaxHeatingWidth( 0.0 );
    _settings->SetPhysicalProblem( problem );
    _settings->SetTimeIntegrator( HEUN );

    if( _settings->GetPhysicalProblem() == EULER1DIDGAS )
        _physicalProblem = new EulerEquationsIdGas( _settings );
    else if( _settings->GetPhysicalProblem() == BURGERS )
        _physicalProblem = new BurgersEquation( _settings );
    else if( _settings->GetPhysicalProblem() == ADVECTION )
        _physicalProblem = new Advection( _settings );

    _mesh = new Mesh1D( _settings->GetXLeft(), _settings->GetXRight(), _settings->GetDimX(), _physicalProblem, _settings );
    _rhs  = new Rhs1D( _settings, _physicalProblem, _mesh->GetDx() );

    Vector* state = new Vector( _physicalProblem->GetStateDim() );
    for( int i = 0; i < _settings->GetDimX(); ++i ) {
        if( _settings->GetPhysicalProblem() == EULER1DIDGAS || _settings->GetPhysicalProblem() == EULER1DREALFLUID ) {
            _physicalProblem->TransformToConservatives( _initialCellStates[i][0], _initialCellStates[i][1], _initialCellStates[i][2], state );
        }
        else {
            for( int j = 0; j < _physicalProblem->GetStateDim(); j++ ) {

                ( *state )[j] = _initialCellStates[i][j];
                // stetige AB zum testen von periodischen RB
                //  double pi     = 3.14159265359;
                //  ( *state )[j] = 2.0 + sin( 2 * pi * ( (double)i * ( _settings->GetXRight() - _settings->GetXLeft() ) /
                //  _settings->GetDimX() ) );
            }
        }
        _mesh->GetCell( i )->SetState( state );
    }

    std::vector<double> conservativeStateMin( _physicalProblem->GetStateDim(), 0 ), conservativeStateMax( _physicalProblem->GetStateDim(), 0 ),
        primitiveStateMin( _physicalProblem->GetStateDim(), 0 ), primitiveStateMax( _physicalProblem->GetStateDim(), 0 );

    for( int j = 0; j < _physicalProblem->GetStateDim(); ++j ) {
        // set min/max with first cell value to adjust range
        conservativeStateMin[j] = ( *( _mesh->GetCell( 0 )->GetState() ) )[j];
        conservativeStateMax[j] = conservativeStateMin[j];
        primitiveStateMin[j]    = _initialCellStates[0][j];
        primitiveStateMax[j]    = primitiveStateMin[j];
    }

    for( int i = 0; i < _settings->GetDimX(); ++i ) {
        for( int j = 0; j < _physicalProblem->GetStateDim(); ++j ) {
            conservativeStateMin[j] = std::min( conservativeStateMin[j], ( *( _mesh->GetCell( i )->GetState() ) )[j] );
            conservativeStateMax[j] = std::max( conservativeStateMax[j], ( *( _mesh->GetCell( i )->GetState() ) )[j] );
            primitiveStateMin[j]    = std::min( primitiveStateMin[j], _initialCellStates[i][j] );
            primitiveStateMax[j]    = std::max( primitiveStateMax[j], _initialCellStates[i][j] );
        }
    }
    if( _settings->GetTimeIntegrator() == HEUN )
        _timeIntegrator = new Heun( _settings, _mesh, _physicalProblem, _rhs );
    else
        exit( EXIT_FAILURE );

    if( _settings->GetTimeIntegrator() == HEUN ) {
        for( unsigned int i = 0; i < _settings->GetNTimesteps(); i++ ) {
            _timeIntegrator->Integrate();
        }
    }

    double* res = new double[( nStates + 1 ) * nCells];

    // store mesh
    for( int i = 0; i < nCells; i++ ) {
        res[i] = _mesh->GetXLeft() + _mesh->GetDx() * ( i + 0.5 );
    }
    // store states
    for( int j = 0; j < nStates; j++ ) {
        for( int i = 0; i < nCells; i++ ) {
            res[i + j * nCells + nCells] = ( *( _mesh->GetCell( i )->GetState() ) )[j];
        }
    }
    plhs[0] = mxCreateNumericMatrix( nCells, nStates + 1, mxDOUBLE_CLASS, mxREAL );
    std::memcpy( mxGetData( plhs[0] ), res, ( nStates + 1 ) * nCells * sizeof( double ) );

    delete res;
    delete _timeIntegrator;
    delete _mesh;
    delete _rhs;
    delete _physicalProblem;
}
