#include "momentsolver.h"
#include <mpi.h>

MomentSolver::MomentSolver( Settings* settings, Mesh* mesh, Problem* problem ) : _settings( settings ), _mesh( mesh ), _problem( problem ) {
    _log         = spdlog::get( "event" );
    _nCells      = _settings->GetNumCells();
    _nMoments    = _settings->GetNMoments();
    _tEnd        = _settings->GetTEnd();
    _nStates     = _settings->GetNStates();
    _nQuadPoints = _settings->GetNQuadPoints();
    _nQTotal     = _settings->GetNQTotal();
    _nTotal      = _settings->GetNTotal();

    _closure = Closure::Create( _settings );
    _time    = TimeSolver::Create( _settings, _mesh );

    _dt = _time->GetTimeStepSize();
}

MomentSolver::~MomentSolver() {
    delete _closure;
    delete _time;
}

void MomentSolver::Solve() {
    auto log = spdlog::get( "event" );

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    // create solution fields
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    if( _settings->HasContinueFile() ) {
        auto ICs = _mesh->Import();
        // TODO: _mesh->Import just imports mean val
    }
    else {
        u = SetupIC();
    }
    MatVec uNew = u;
    MatVec uOld = u;
    MatVec uQ   = MatVec( _nCells + 1, Matrix( _nStates, _settings->GetNqPE() ) );
    _lambda     = MatVec( _nCells + 1, Matrix( _nStates, _nTotal ) );

    Vector ds( _nStates );
    Vector u0( _nStates );

    double residualFull;

    std::vector<unsigned> cellIndexPE = _settings->GetCellIndexPE();
    std::vector<int> PEforCell        = _settings->GetPEforCell();

    log->info( "PE {0}: kStart {1}, kEnd {2}", _settings->GetMyPE(), _settings->GetKStart(), _settings->GetKEnd() );

    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            u0[l] = u[j]( l, 0 );
        }
        _closure->DS( ds, u0 );
        for( unsigned l = 0; l < _nStates; ++l ) {
            _lambda[j]( l, 0 ) = ds[l];
        }
    }

    // Converge initial condition entropy variables for One Shot IPM
    if( _settings->GetMaxIterations() == 1 ) {
        _settings->SetMaxIterations( 1000 );
        for( unsigned j = 0; j < _nCells; ++j ) _closure->SolveClosure( _lambda[j], u[j] );
        _settings->SetMaxIterations( 1 );
    }

    auto numFluxPtr = std::bind( &MomentSolver::numFlux,
                                 this,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5 );

    log->info( "{:10}   {:10}", "t", "residual" );
    // Begin time loop
    for( double t = 0.0; t < _tEnd; t += _dt ) {
        double residual = 0;

#pragma omp parallel for schedule( dynamic, 10 )
        for( unsigned j = 0; j < cellIndexPE.size(); ++j ) {
            _closure->SolveClosure( _lambda[cellIndexPE[j]], u[cellIndexPE[j]] );
        }

        // MPI Broadcast lambdas to all PEs
        for( unsigned j = 0; j < _nCells; ++j ) {
            uOld[j] = u[j];    // save old Moments for residual computation
            MPI_Bcast( _lambda[j].GetPointer(), int( _nStates * _nTotal ), MPI_DOUBLE, PEforCell[j], MPI_COMM_WORLD );
        }

        for( unsigned j = 0; j < _nCells; ++j ) {
            uQ[j] = _closure->U( _closure->EvaluateLambdaOnPE( _lambda[j] ) );
            u[j].reset();
            multOnPENoReset( uQ[j], _closure->GetPhiTildeWf(), u[j], _settings->GetKStart(), _settings->GetKEnd() );
        }

        _time->Advance( numFluxPtr, uNew, u, uQ );

        // perform reduction to obtain full moments
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Reduce( uNew[j].GetPointer(), u[j].GetPointer(), int( _nStates * _nTotal ), MPI_DOUBLE, MPI_SUM, PEforCell[j], MPI_COMM_WORLD );
        }
        for( unsigned j = 0; j < cellIndexPE.size(); ++j ) {
            residual += std::fabs( u[cellIndexPE[j]]( 0, 0 ) - uOld[cellIndexPE[j]]( 0, 0 ) ) * _mesh->GetArea( cellIndexPE[j] ) / _dt;
        }
        MPI_Reduce( &residual, &residualFull, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        if( _settings->GetMyPE() == 0 ) log->info( "{:03.8f}   {:01.5e}", t, residualFull );
    }
    if( _settings->GetMyPE() != 0 ) return;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    log->info( "" );
    log->info( "Finished!" );
    log->info( "" );
    log->info( "Runtime: {0}s", std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 );

    Matrix meanAndVar( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    Matrix phiTildeWf = _closure->GetPhiTildeWf();
    Vector tmp( _nStates, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i, j ) += tmp[i] * phiTildeWf( k, 0 );
            }
        }

        // var
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i + _nStates, j ) += pow( tmp[i] - meanAndVar( i, j ), 2 ) * phiTildeWf( k, 0 );
            }
        }
    }

    if( _settings->GetProblemType() == P_SHALLOWWATER_2D )
        _mesh->ExportShallowWater( meanAndVar );
    else
        _mesh->Export( meanAndVar );

    unsigned evalCell  = 0;    // 2404;
    unsigned plotState = 0;
    _mesh->PlotInXi( _closure->U( _closure->EvaluateLambda( _lambda[evalCell] ) ), plotState );
}

void MomentSolver::numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    // out += _problem->G( u1, u2, nUnit, n ) * _closure->GetPhiTildeWf();
    multOnPENoReset( _problem->G( u1, u2, nUnit, n ), _closure->GetPhiTildeWf(), out, _settings->GetKStart(), _settings->GetKEnd() );
}

void MomentSolver::CalculateMoments( MatVec& out, const MatVec& lambda ) {
    Matrix U( _nStates, _nQTotal, 0.0 );
    Matrix evalLambda( _nStates, _nQTotal, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        _closure->EvaluateLambda( evalLambda, lambda[j] );
        _closure->U( U, evalLambda );
        out[j] = U * _closure->GetPhiTildeWf();
    }
}

MatVec MomentSolver::SetupIC() {
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    std::vector<Polynomial*> quad = _closure->GetQuadrature();
    Vector xiEta( _settings->GetNDimXi() );
    Matrix uIC( _nStates, _nQTotal );
    Matrix phiTildeWf = _closure->GetPhiTildeWf();
    unsigned n;
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index =
                    unsigned( ( k - k % unsigned( std::pow( _nQuadPoints, l ) ) ) / unsigned( std::pow( _nQuadPoints, l ) ) ) % _nQuadPoints;
                xiEta[l] = quad[n]->GetNodes()[index];
            }
            // std::cout << xiEta << std::endl;
            column( uIC, k ) = IC( _mesh->GetCenterPos( j ), xiEta );
        }
        u[j] = uIC * phiTildeWf;
    }
    return u;
}

Vector MomentSolver::IC( Vector x, Vector xi ) {
    Vector y( _nStates );
    if( _settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
        if( xi.size() == 1 ) {
            double a     = 0.5;
            double b     = 1.5;
            double sigma = 0.2;    // 0.2
            double uL    = 12.0;
            double uR    = 3.0;
            if( x[0] < a + sigma * xi[0] ) {
                y[0] = uL;
                return y;
            }
            else if( x[0] < b + sigma * xi[0] ) {
                y[0] = uL + ( uR - uL ) * ( a + sigma * xi[0] - x[0] ) / ( a - b );
                return y;
            }
            else {
                y[0] = uR;
                return y;
            }
        }
        else if( xi.size() == 2 ) {
            double x0     = 0.3;
            double x1     = 0.6;
            double sigma0 = 0.2;    // 0.2
            double sigma1 = 0.1;
            double uL     = 12.0;
            double uM     = 6.0;
            double uR     = 1.0;

            if( x[0] < x0 )
                y[0] = uL + sigma0 * xi[0];
            else if( x[0] < x1 )
                y[0] = uM + sigma1 * xi[1];
            else
                y[0] = uR;
            return y;
        }
    }
    if( _settings->GetProblemType() == ProblemType::P_SHALLOWWATER_1D ) {
        if( xi.size() == 1 ) {
            double a     = 1000.0;
            double sigma = 0.0;    // 0.2
            double uL    = 10.0;
            double hLhR  = 0.05;
            double uR    = uL * hLhR;
            y[1]         = 0.0;
            if( x[0] < a + sigma * xi[0] ) {
                y[0] = uL;
                return y;
            }
            else {
                y[0] = uR;
                return y;
            }
        }
        else if( xi.size() == 2 ) {
            double x0     = 0.3;
            double x1     = 0.6;
            double sigma0 = 0.2;    // 0.2
            double sigma1 = 0.1;
            double uL     = 12.0;
            double uM     = 6.0;
            double uR     = 1.0;

            if( x[0] < x0 )
                y[0] = uL + sigma0 * xi[0];
            else if( x[0] < x1 )
                y[0] = uM + sigma1 * xi[1];
            else
                y[0] = uR;
            return y;
        }
    }
    if( _settings->GetProblemType() == ProblemType::P_SHALLOWWATER_2D ) {
        double a      = 0.5;
        double sigma  = 1.0;    // 0.5
        double sigma1 = 0.5;
        double uL     = 10.0;    // 10.0;
        double uR     = 5.0;
        y[1]          = 0.0;
        y[2]          = 0.0;
        if( xi.size() == 1 ) {
            if( x[0] < a ) {
                y[0] = uL + sigma * xi[0];
                return y;
            }
            else {
                y[0] = uR;    // - sigma * xi[0];
                return y;
            }
        }
        else if( xi.size() == 2 ) {
            if( x[0] < a ) {
                y[0] = uL + sigma * xi[0];
                return y;
            }
            else {
                y[0] = uR + sigma1 * xi[1];
                return y;
            }
        }
    }
    else if( _settings->GetProblemType() == ProblemType::P_EULER_1D ) {
        double x0    = 0.3;
        double sigma = 0.0;
        double gamma = 1.4;

        double rhoL = 1.0;
        double rhoR = 0.125;
        double pL   = 1.0;
        double pR   = 0.1;
        double uL   = 0.0;
        double uR   = 0.0;
        if( x[0] < x0 + sigma * xi[0] ) {
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
    else if( _settings->GetProblemType() == ProblemType::P_EULER_2D && false ) {    // pipe testcase
        double sigma = 0.2;
        double gamma = 1.4;
        double R     = 287.87;
        double T     = 273.15;
        double p     = 101325.0;
        double Ma    = 0.0;
        double a     = sqrt( gamma * R * T );
        double pi    = 3.14159265359;

        double uMax  = Ma * a;
        double angle = ( 1.25 + sigma * 0.0 ) * ( 2.0 * pi ) / 360.0;
        double uF    = uMax * cos( angle );
        double vF    = uMax * sin( angle );

        double rhoFarfield = p / ( R * T );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        /*
        if( x[1] < 1.1 + sigma * xi[0]  ) {
            y[0] = 0.5 * rhoFarfield;
            y[1] = rhoFarfield * uF;
            y[2] = rhoFarfield * vF;
            y[3] = 0.5 * ( kineticEnergyL + innerEnergyL );
        }*/
        if( x[1] < 1.1 + sigma && x[1] < 1.1 - sigma ) {
            y[0] = 0;
        }
        if( x[1] < 1.1 + sigma && x[1] > 1.1 - sigma ) {
            y[0] = 1;
        }
        if( x[1] > 1.1 + sigma && x[1] > 1.1 - sigma ) {
            y[0] = 2;
        }
        return y;
    }
    else if( _settings->GetProblemType() == ProblemType::P_EULER_2D && true ) {
        double sigma  = 0.5;
        double sigma1 = 0.01;
        double gamma  = 1.4;
        double R      = 287.87;
        double T      = 273.15;
        double p      = 101325.0;
        double Ma     = 0.8;
        if( xi.size() == 2 ) {
            Ma = Ma + xi[1] * sigma1;
        }
        double a  = sqrt( gamma * R * T );
        double pi = 3.14159265359;

        double uMax  = Ma * a;
        double angle = ( 1.25 + sigma * xi[0] ) * ( 2.0 * pi ) / 360.0;
        double uF    = uMax * cos( angle );
        double vF    = uMax * sin( angle );

        double rhoFarfield = p / ( R * T );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        return y;
    }
    _log->error( "Reached end of IC. No initial condition set" );
    exit( EXIT_FAILURE );
}
