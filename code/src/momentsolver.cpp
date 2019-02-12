#include "momentsolver.h"
#include <mpi.h>

MomentSolver::MomentSolver( Settings* settings, Mesh* mesh, Problem* problem ) : _settings( settings ), _mesh( mesh ), _problem( problem ) {
    _log         = spdlog::get( "event" );
    _nCells      = _settings->GetNumCells();
    _nMoments    = _settings->GetNMoments();
    _tStart      = 0.0;
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
    if( _settings->HasRestartFile() ) {
        u = this->Import();
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
        // TODO: Recalculate moments here
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
    double t;
    for( t = _tStart; t < _tEnd; t += _dt ) {
        double residual = 0;

#pragma omp parallel for schedule( dynamic, 10 )
        for( unsigned j = 0; j < static_cast<unsigned>( cellIndexPE.size() ); ++j ) {
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

    Matrix meanAndVar;
    if( _settings->HasExactSolution() ) {
        meanAndVar = Matrix( 4 * _nStates, _mesh->GetNumCells(), 0.0 );
    }
    else {
        meanAndVar = Matrix( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    }
    Matrix phiTildeWf = _closure->GetPhiTildeWf();
    Vector tmp( _nStates, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        // expected value
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i, j ) += tmp[i] * phiTildeWf( k, 0 );
            }
        }

        // variance
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i + _nStates, j ) += pow( tmp[i] - meanAndVar( i, j ), 2 ) * phiTildeWf( k, 0 );
            }
        }
        if( _settings->HasExactSolution() ) {
            Vector xiEta( _settings->GetNDimXi() );
            std::vector<Polynomial*> quad = _closure->GetQuadrature();
            unsigned n;

            for( unsigned k = 0; k < _nQTotal; ++k ) {
                for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                    if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                    if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                    unsigned index =
                        unsigned( ( k - k % unsigned( std::pow( _nQuadPoints, l ) ) ) / unsigned( std::pow( _nQuadPoints, l ) ) ) % _nQuadPoints;
                    xiEta[l] = quad[n]->GetNodes()[index];
                }
                tmp = _problem->ExactSolution( t, _mesh->GetCenterPos( j ), xiEta );
                // expected value exact
                for( unsigned i = 0; i < _nStates; ++i ) {
                    meanAndVar( 2 * _nStates + i, j ) += tmp[i] * phiTildeWf( k, 0 );
                }
                // variance exact
                for( unsigned i = 0; i < _nStates; ++i ) {
                    meanAndVar( 3 * _nStates + i, j ) += pow( tmp[i] - meanAndVar( 2 * _nStates + i, j ), 2 ) * phiTildeWf( k, 0 );
                }
            }
        }
    }

    if( _settings->GetProblemType() == P_SHALLOWWATER_2D )
        _mesh->ExportShallowWater( meanAndVar );
    else
        _mesh->Export( meanAndVar );

    this->Export( u );

    unsigned evalCell  = 300;    // 2404;
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
    std::vector<Vector> IC;
    if( _settings->HasICFile() ) {
        IC = _mesh->Import();
    }
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

            if( _settings->HasICFile() ) {
                column( uIC, k ) = _problem->LoadIC( IC[j], xiEta );
            }
            else {
                column( uIC, k ) = _problem->IC( _mesh->GetCenterPos( j ), xiEta );
            }
        }

        u[j] = uIC * phiTildeWf;
    }
    return u;
}

void MomentSolver::Export( const MatVec& u ) const {
    std::shared_ptr<spdlog::logger> writer = spdlog::get( "moments" );
    writer->info( "{0}", _tEnd );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned j = 0; j < _nStates; ++j ) {

            for( unsigned k = 0; k < _nTotal - 1; ++k ) {
                line << u[i]( j, k ) << ",";
            }
        }
        line << u[i]( _nStates - 1, _nTotal - 1 );
        writer->info( line.str() );
    }
    writer->flush();
}

MatVec MomentSolver::Import() {
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() );
    std::string line;
    std::getline( file, line );
    _tStart = std::stod( line );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < _nTotal; ++k ) {
                std::getline( lineStream, cell, ',' );
                u[i]( j, k ) = std::stod( cell );
            }
        }
    }
    file.close();
    return u;
}
