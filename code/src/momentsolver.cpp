﻿#include "momentsolver.h"

MomentSolver::MomentSolver( Problem* problem ) : _problem( problem ) {
    _quad = new Legendre( _problem->GetNQuadPoints() );

    _mesh   = _problem->GetMesh();
    _nCells = _mesh->GetNumCells();
    _a      = 0.0;
    _b      = 3.0;

    _uL = 12.0;
    _uR = 3.0;

    _nMoments = _problem->GetNMoments();
    _tEnd     = _problem->GetTEnd();
    _nStates  = _problem->GetNStates();

    _closure = Closure::Create( _problem );
    _limiter = Limiter::Create( _closure, _problem );
    _time    = TimeSolver::Create( _problem, _closure );

    _dt = _time->GetTimeStepSize();

    _plotEngine = new PlotEngine( _problem );
}

MomentSolver::~MomentSolver() {
    delete _quad;
    delete _closure;
    delete _limiter;
    delete _time;
    delete _plotEngine;
}

void MomentSolver::Solve() {
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    double t = 0;
    // create solution fields
    std::vector<Matrix> uNew( _nCells, Matrix( _nStates, _nMoments, 0.0 ) );
    _lambda                = std::vector<Matrix>( _nCells + 1, Matrix( _nStates, _nMoments, 0.0 ) );
    std::vector<Matrix> u  = SetupIC();
    std::vector<Matrix> uQ = std::vector<Matrix>( _nCells + 1, Matrix( _nStates, _problem->GetNQuadPoints(), 0.0 ) );

    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _problem->GetProblemType() == "Euler" ) {
            _lambda[j]( 2, 0 ) = -1.0;
            _lambda[j]( 0, 0 ) = 1.0;
        }
        if( _problem->GetProblemType() == "Euler2D" ) {
            _lambda[j]( 3, 0 ) = -1.0;
            _lambda[j]( 0, 0 ) = 1.0;
        }
        _lambda[j] = _closure->SolveClosure( u[j], _lambda[j] );
        u[j]       = CalculateMoments( _lambda[j] );    // kann raus!
    }

    // Begin time loop
    while( t < _tEnd ) {
        // Modify moments into realizable direction
        for( unsigned j = 0; j < _nCells; ++j ) {
            u[j]  = CalculateMoments( _lambda[j] );
            uQ[j] = _closure->U( _closure->EvaluateLambda( _lambda[j] ) );
        }

        // Time Update Moments
        _time->Advance(
            std::bind( &MomentSolver::numFlux, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 ),
            uNew,
            u,
            uQ );

        // Time Update dual variables
        //#pragma omp parallel for
        for( unsigned j = 0; j < _nCells; ++j ) {
            _lambda[j] = _closure->SolveClosure( uNew[j], _lambda[j] );
        }

        // Plot( t );

        if( std::fabs( t ) <= std::numeric_limits<double>::epsilon() )
            std::cout << std::fixed << std::setprecision( 8 ) << "t = " << t << std::flush;
        else
            std::cout << std::fixed << std::setprecision( 8 ) << "\r"
                      << "t = " << t << std::flush;
        t += _dt;
    }

    _tEnd = t;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    std::cout << "\nFinished!\nRuntime: " << std::setprecision( 3 )
              << std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 << "s" << std::endl;

    Matrix res( _nStates, _mesh->GetNumCells() );
    Vector tmp( _nStates, 0.0 );
    auto xiQuad = _quad->GetNodes();
    Vector w    = _quad->GetWeights();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned i = 0; i < _nStates; ++i ) {
            for( unsigned k = 0; k < _problem->GetNQuadPoints(); ++k ) {
                _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                res( i, j ) += 0.5 * w[k] * tmp[i];
            }
        }
    }
    _mesh->Export( res );
}

Matrix MomentSolver::numFlux( const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    Matrix g = _problem->G( u1, u2, nUnit, n );
    return 0.5 * g * _closure->GetPhiTildeW();
}

Matrix MomentSolver::CalculateMoments( const Matrix& lambda ) {
    return 0.5 * _closure->U( _closure->EvaluateLambda( lambda ) ) * _closure->GetPhiTildeW();
}

std::vector<Matrix> MomentSolver::SetupIC() {
    std::vector<Matrix> out( _nCells + _problem->GetMesh()->GetNBoundaries(), Matrix( _nStates, _nMoments, 0.0 ) );
    Vector xi = _quad->GetNodes();
    Matrix uIC( _nStates, _problem->GetNQuadPoints(), 0.0 );
    Matrix phiTildeW = _closure->GetPhiTildeW();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _problem->GetNQuadPoints(); ++k ) {
            column( uIC, k ) = IC( _problem->GetMesh()->GetCenterPos( j ), xi[k] );
        }
        out[j] = 0.5 * uIC * phiTildeW;
    }
    return out;
}

Vector MomentSolver::IC( Vector x, double xi ) {
    Vector y( _nStates );
    if( _problem->GetProblemType() == "Burgers" ) {
        double a     = 0.5;
        double b     = 1.5;
        double sigma = 0.2;
        if( x[0] < a + sigma * xi ) {
            y[0] = _uL;
            return y;
        }
        else if( x[0] < b + sigma * xi ) {
            y[0] = _uL + ( _uR - _uL ) * ( a + sigma * xi - x[0] ) / ( a - b );
            return y;
        }
        else {
            y[0] = _uR;
            return y;
        }
    }
    else if( _problem->GetProblemType() == "Euler" ) {
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
    else if( _problem->GetProblemType() == "Euler2D" ) {
        double sigma = 0.0;
        double gamma = 1.4;

        double rhoFarfield = 1.0;
        double pFarfield   = 1.0;
        double uMax        = 1.0;
        double angle       = 0.0 + sigma;
        double uF          = uMax * cos( angle );
        double vF          = uMax * sin( angle );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( pFarfield / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        if( x[0] > 0.0 ) {
            rhoFarfield = 0.8;
            pFarfield   = 0.8;
            uMax        = 1.0;
            angle       = 0.0 + sigma;
            uF          = uMax * cos( angle );
            vF          = uMax * sin( angle );
            // y[0]           = rhoFarfield;
            // y[1]           = rhoFarfield * uF;
            // y[2]           = rhoFarfield * vF;
            kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
            innerEnergyL   = ( pFarfield / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
            // y[3]           = kineticEnergyL + innerEnergyL;
        }
        return y;
    }
    std::cerr << "Reached end of IC. No initial condition set" << std::endl;
    exit( EXIT_FAILURE );
}

void MomentSolver::Plot( double time ) {

    static unsigned plotCtr = 0;
    if( _problem->GetMesh()->GetDimension() == 1 ) {
        const unsigned int nQuadFine = 200;
        Legendre quadFine( nQuadFine );
        Vector wFine      = quadFine.GetWeights();
        Vector xiQuadFine = quadFine.GetNodes();
        Vector w          = _quad->GetWeights();
        Vector xiQuad     = _quad->GetNodes();
        Vector x          = _mesh->GetNodePositionsX();
        unsigned nFine;

        try {
            auto file = cpptoml::parse_file( _problem->GetInputFile() );
            auto plot = file->get_table( "plot" );

            nFine = plot->get_as<unsigned>( "evalPoints" ).value_or( 0 );
        } catch( const cpptoml::parse_exception& e ) {
            std::cerr << "Failed to parse " << _problem->GetInputFile() << ": " << e.what() << std::endl;
            exit( EXIT_FAILURE );
        }

        Vector xFine( nFine, 0.0 ), exResMean( nFine, 0.0 ), resMean( _nCells, 0.0 );

        for( unsigned int j = 0; j < nFine; ++j ) {
            xFine[j] = _a + j * ( _b - _a ) / ( nFine - 1 );
            for( unsigned int k = 0; k < nQuadFine; ++k ) {
                exResMean[j] += 0.5 * wFine[k] * _problem->ExactSolution( _tEnd, xFine[j], xiQuadFine[k] );
            }
        }
        Vector resVec( _nStates, 0.0 );
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned k = 0; k < _problem->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                resMean[j] += 0.5 * w[k] * resVec[0];
            }
        }

        Vector exResVar( nFine, 0.0 ), resVar( _nCells, 0.0 );
        for( unsigned int j = 0; j < nFine; ++j ) {
            for( unsigned int k = 0; k < nQuadFine; ++k ) {
                exResVar[j] += 0.5 * wFine[k] * _problem->ExactSolution( _tEnd, xFine[j], xiQuadFine[k] );
            }
        }

        resVec.reset();
        double expectValue;
        for( unsigned j = 0; j < _nCells; ++j ) {
            expectValue = 0.0;
            for( unsigned k = 0; k < _problem->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                expectValue += 0.5 * w[k] * resVec[0];
            }
            for( unsigned k = 0; k < _problem->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
                resVar[j] += 0.5 * w[k] * pow( resVec[0] - expectValue, 2 );
            }
        }

        if( plotCtr == 0 ) {
            _plotEngine->setupPlot( 0, "Mean", "x", "y" );
            _plotEngine->setupPlot( 1, "Var", "x", "y" );
            _plotEngine->addPlotData( 0, Result1D{x, resMean}, "simulation" );
            //_plotEngine->addPlotData( 0, Result1D{xFine, exResMean}, "exact" );
            _plotEngine->addPlotData( 1, Result1D{x, resVar}, "simulation" );
            //_plotEngine->addPlotData( 1, Result1D{xFine, exResVar}, "exact" );
        }
        else {
            _plotEngine->updatePlotData( 0, Result1D{x, resMean}, "simulation" );
            //            _plotEngine->updatePlotData( 0, Result1D{xFine, exResMean}, "exact" );
            _plotEngine->updatePlotData( 1, Result1D{x, resVar}, "simulation" );
            //_plotEngine->updatePlotData( 1, Result1D{_x, exResVar}, "exact" );
        }
        _plotEngine->replot();
        plotCtr++;
    }
    else if( _problem->GetMesh()->GetDimension() == 2 ) {
        // TODO
    }
    if( plotCtr == 1 ) {
        _plotEngine->show();
    }
}
