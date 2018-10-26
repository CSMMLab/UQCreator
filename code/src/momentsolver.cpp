﻿#include "momentsolver.h"

MomentSolver::MomentSolver( Settings* settings, Mesh* mesh, Problem* problem ) : _settings( settings ), _mesh( mesh ), _problem( problem ) {
    _nCells      = _settings->GetNumCells();
    _nMoments    = _settings->GetNMoments();
    _tEnd        = _settings->GetTEnd();
    _nStates     = _settings->GetNStates();
    _nQuadPoints = _settings->GetNQuadPoints();

    _quad    = new Legendre( _nQuadPoints );
    _closure = Closure::Create( _settings );
    _time    = TimeSolver::Create( _settings, _mesh );

    _dt = _time->GetTimeStepSize();

    if( _settings->GetPlotEnabled() ) {
        _plotEngine = new PlotEngine( _settings, _closure, _mesh, _problem );
    }
    else {
        _plotEngine = nullptr;
    }
}

MomentSolver::~MomentSolver() {
    delete _quad;
    delete _closure;
    delete _time;
    delete _plotEngine;
}

void MomentSolver::Solve() {
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    double t = 0;
    // create solution fields
    MatVec uNew( _nCells, Matrix( _nStates, _nMoments, 0.0 ) );
    _lambda = MatVec( _nCells + 1, Matrix( _nStates, _nMoments, 0.0 ) );
    MatVec u( _nCells );
    SetupIC( u );
    MatVec uQ = MatVec( _nCells + 1, Matrix( _nStates, _nQuadPoints, 0.0 ) );

    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _settings->GetProblemType() == ProblemType::P_EULER_1D || _settings->GetProblemType() == ProblemType::P_EULER_2D ) {
            _lambda[j]( 0, 0 )                           = 1.0;
            _lambda[j]( _settings->GetNStates() - 1, 0 ) = -1.0;
        }
        _closure->SolveClosure( _lambda[j], u[j] );
    }

    // Begin time loop
    unsigned nSteps = 0;
    while( t < _tEnd ) {
        // Modify moments into realizable direction
        CalculateMoments( u, _lambda );
        for( unsigned j = 0; j < _nCells; ++j ) {
            _closure->U( uQ[j], _closure->EvaluateLambda( _lambda[j] ) );
        }

        // Time Update Moments
        _time->Advance( std::bind( &MomentSolver::numFlux,
                                   this,
                                   std::placeholders::_1,
                                   std::placeholders::_2,
                                   std::placeholders::_3,
                                   std::placeholders::_4,
                                   std::placeholders::_5 ),
                        uNew,
                        u,
                        uQ );

        // Time Update dual variables
        //#pragma omp parallel for
        for( unsigned j = 0; j < _nCells; ++j ) {
            _closure->SolveClosure( _lambda[j], uNew[j] );
        }

        Plot( t, nSteps );
        double residual = 0;
        for( unsigned j = 0; j < _nCells; ++j ) {
            residual += std::fabs( uNew[j]( 0, 0 ) - u[j]( 0, 0 ) );
        }

        std::cout << std::fixed << std::setprecision( 8 ) << "\r"
                  << "t = " << t << " -> E[rho] residual is " << residual << std::endl;
        t += _dt;
        nSteps++;
    }

    _tEnd = t;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    std::cout << "\nFinished!\nRuntime: " << std::setprecision( 3 )
              << std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 << "s" << std::endl;

    Matrix meanAndVar( 2 * _nStates, _mesh->GetNumCells() );
    Vector tmp( _nStates, 0.0 );
    auto xiQuad = _quad->GetNodes();
    Vector w    = _quad->GetWeights();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned i = 0; i < _nStates; ++i ) {
            // mean
            for( unsigned k = 0; k < _nQuadPoints; ++k ) {
                _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                meanAndVar( i, j ) += 0.5 * w[k] * tmp[i];
            }
            // var
            for( unsigned k = 0; k < _nQuadPoints; ++k ) {
                _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                meanAndVar( i + _nStates, j ) += 0.5 * w[k] * pow( tmp[i] - meanAndVar( i, j ), 2 );
            }
        }
    }

    _mesh->Export( meanAndVar );
}

void MomentSolver::numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    out += 0.5 * _problem->G( u1, u2, nUnit, n ) * _closure->GetPhiTildeW();
}

void MomentSolver::CalculateMoments( MatVec& out, const MatVec& lambda ) {
    Matrix U( _nStates, _nQuadPoints, 0.0 );
    Matrix evalLambda( _nStates, _nQuadPoints, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        _closure->EvaluateLambda( evalLambda, lambda[j] );
        _closure->U( U, evalLambda );
        out[j] = 0.5 * U * _closure->GetPhiTildeW();
    }
}

void MomentSolver::SetupIC( MatVec& out ) {
    Vector xi = _quad->GetNodes();
    Matrix uIC( _nStates, _nQuadPoints, 0.0 );
    Matrix phiTildeW = _closure->GetPhiTildeW();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            column( uIC, k ) = IC( _mesh->GetCenterPos( j ), xi[k] );
        }
        out[j] = 0.5 * uIC * phiTildeW;
    }
}

Vector MomentSolver::IC( Vector x, double xi ) {
    Vector y( _nStates );
    if( _settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
        double a     = 0.5;
        double b     = 1.5;
        double sigma = 0.2;
        double uL    = 12.0;
        double uR    = 3.0;
        if( x[0] < a + sigma * xi ) {
            y[0] = uL;
            return y;
        }
        else if( x[0] < b + sigma * xi ) {
            y[0] = uL + ( uR - uL ) * ( a + sigma * xi - x[0] ) / ( a - b );
            return y;
        }
        else {
            y[0] = uR;
            return y;
        }
    }
    else if( _settings->GetProblemType() == ProblemType::P_EULER_1D ) {
        double x0    = 0.3;
        double sigma = 0.05;
        double gamma = 1.4;

        double rhoL = 1.0;
        double rhoR = 0.3;
        double pL   = 1.0;
        double pR   = 0.3;
        double uL   = 0.0;
        double uR   = 0.0;
        if( x[0] < x0 + sigma * xi ) {
            y[0]                  = rhoL;
            y[1]                  = rhoL * uL;
            double kineticEnergyL = 0.5 * rhoL * pow( uL, 2 );
            double innerEnergyL   = ( pL / ( rhoL * ( gamma - 1 ) ) ) * rhoL;
            y[2]                  = kineticEnergyL + innerEnergyL;
        }
        else {
            y[0]                  = rhoR;
            y[1]                  = rhoR * uR;
            double kineticEnergyR = 0.5 * rhoR * pow( uR, 2 );
            double innerEnergyR   = ( pR / ( rhoR * ( gamma - 1 ) ) ) * rhoR;
            y[2]                  = kineticEnergyR + innerEnergyR;
        }
        return y;
    }
    else if( _settings->GetProblemType() == ProblemType::P_EULER_2D ) {
        double sigma = 0.1;
        double gamma = 1.4;

        double rhoFarfield = 1.0;
        double pFarfield   = 1.0;
        double uMax        = 0.1;
        double angle       = 0.3 + sigma * xi;
        double uF          = uMax * cos( angle );
        double vF          = uMax * sin( angle );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( pFarfield / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        return y;
    }
    std::cerr << "Reached end of IC. No initial condition set" << std::endl;
    exit( EXIT_FAILURE );
}

void MomentSolver::Plot( double time, unsigned nSteps ) {
    static double lastPlotTime = 0.0;
    if( _plotEngine != nullptr ) {
        bool execPlot = false;
        if( std::fabs( _settings->GetPlotTimeInterval() + 1.0 ) < std::numeric_limits<double>::epsilon() &&
            time - lastPlotTime > _settings->GetPlotTimeInterval() ) {
            lastPlotTime = time;
            execPlot     = true;
        }
        if( _settings->GetPlotStepInterval() != 0 && nSteps % _settings->GetPlotStepInterval() == 0 ) {
            execPlot = true;
        }
        if( execPlot ) {
            _plotEngine->updatePlotData( time, _lambda );
        }
        else {
            _plotEngine->refresh();
        }
    }
}
