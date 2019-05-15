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
    _settings->SetDT( _dt );

    if( _settings->HasReferenceFile() ) {
        _referenceSolution = _mesh->Import();
        assert( _nCells == _referenceSolution.size() );
    }
}

MomentSolver::~MomentSolver() {
    delete _closure;
    delete _time;
}

void MomentSolver::Solve() {
    bool useAdaptivity = true;                                                    // flag for using adaptivity
    VectorU refinementLevel( _nCells, _settings->GetNRefinementLevels() - 1 );    // vector carries refinement level for each cell
    Matrix refinementIndicatorPlot( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    VectorU nTotal = _settings->GetNTotalRefinementLevel();

    auto log                                  = spdlog::get( "event" );
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    // create settings+closure for restart as well as solution fields
    Settings* prevSettings = DeterminePreviousSettings();
    Closure* prevClosure   = DeterminePreviousClosure( prevSettings );
    if( _settings->HasRestartFile() ) _tStart = prevSettings->GetTEnd();
    MatVec u = DetermineMoments( prevSettings->GetNTotal() );
    SetDuals( prevSettings, prevClosure, u );
    MatVec uQ = MatVec( _nCells + 1, Matrix( _nStates, _settings->GetNqPE() ) );

    std::vector<unsigned> cellIndexPE = _settings->GetCellIndexPE();
    std::vector<int> PEforCell        = _settings->GetPEforCell();

    // log->info( "PE {0}: kStart {1}, kEnd {2}", _settings->GetMyPE(), _settings->GetKStart(), _settings->GetKEnd() );

    MatVec uNew = u;
    MatVec uOld = u;

    // set up function pointer for right hand side
    auto numFluxPtr = std::bind( &MomentSolver::numFlux,
                                 this,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5,
                                 std::placeholders::_6 );

    if( _settings->GetMyPE() == 0 ) log->info( "{:10}   {:10}", "t", "residual" );

    // init time and residual
    double t      = _tStart;
    int timeIndex = 0;
    double dt;
    double minResidual  = _settings->GetMinResidual();
    double residualFull = minResidual + 1.0;

    // Begin time loop
    while( t < _tEnd && residualFull > minResidual ) {
        double residual = 0;
#pragma omp parallel for schedule( dynamic, 10 )
        for( unsigned j = 0; j < static_cast<unsigned>( cellIndexPE.size() ); ++j ) {
            _closure->SolveClosure( _lambda[cellIndexPE[j]], u[cellIndexPE[j]], nTotal[refinementLevel[cellIndexPE[j]]], _nQTotal );
        }

        // std::cout << "Before MPI_Bcast" << std::endl;
        // MPI Broadcast lambdas to all PEs
        for( unsigned j = 0; j < _nCells; ++j ) {
            uOld[j] = u[j];    // save old Moments for residual computation
            MPI_Bcast( _lambda[j].GetPointer(), int( _nStates * _nTotal ), MPI_DOUBLE, PEforCell[j], MPI_COMM_WORLD );
        }

        // std::cout << "Before uQ Computation" << std::endl;
        // for nQ refinement: here we need new refinement level for nQTotal and old level for nTotal
        // compute solution at quad points
        for( unsigned j = 0; j < _nCells; ++j ) {
            uQ[j] = _closure->U( _closure->EvaluateLambdaOnPE( _lambda[j], nTotal[refinementLevel[j]] ) );
        }

        // std::cout << "Before refinement Computation" << std::endl;
        // determine refinement level of cells on current PE
        for( unsigned j = 0; j < static_cast<unsigned>( cellIndexPE.size() ); ++j ) {
            double indicator = std::fabs( u[cellIndexPE[j]]( 0, nTotal[refinementLevel[cellIndexPE[j]]] - 1 ) ) +
                               std::fabs( u[cellIndexPE[j]]( 0, nTotal[refinementLevel[cellIndexPE[j]]] - 2 ) );
            if( indicator > 0.001 && refinementLevel[cellIndexPE[j]] < _settings->GetNRefinementLevels() - 1 )
                refinementLevel[cellIndexPE[j]] += 1;
            else if( indicator < 0.0001 && refinementLevel[cellIndexPE[j]] > 0 )
                refinementLevel[cellIndexPE[j]] -= 1;
        }

        // std::cout << "Before MPI_Bcast" << std::endl;
        // broadcast refinemt level to all PEs
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Bcast( &refinementLevel[j], 1, MPI_UNSIGNED, PEforCell[j], MPI_COMM_WORLD );
        }

        // compute partial moment vectors on each PE (for inexact dual variables)
        for( unsigned j = 0; j < _nCells; ++j ) {
            u[j].reset();
            multOnPENoReset( uQ[j], _closure->GetPhiTildeWf(), u[j], _settings->GetKStart(), _settings->GetKEnd(), nTotal[refinementLevel[j]] );
        }

        // determine time step size
        double dtMinOnPE = 1e10;
        double dtCurrent;
        for( unsigned j = 0; j < _nCells; ++j ) {
            dtCurrent = _problem->ComputeDt( uQ[j], _mesh->GetMaxEdge( j ) / _mesh->GetArea( j ) );
            if( dtCurrent < dtMinOnPE ) dtMinOnPE = dtCurrent;
        }
        MPI_Allreduce( &dtMinOnPE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
        // std::cout << "dt = " << dt << std::endl;
        t += dt;
        ++timeIndex;

        _time->Advance( numFluxPtr, uNew, u, uQ, dt, refinementLevel );

        // perform reduction onto u to obtain full moments on PE PEforCell[j], which is the PE that solves the dual problem for cell j
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Reduce( uNew[j].GetPointer(), u[j].GetPointer(), int( _nStates * _nTotal ), MPI_DOUBLE, MPI_SUM, PEforCell[j], MPI_COMM_WORLD );
        }

        // compute residual
        for( unsigned j = 0; j < cellIndexPE.size(); ++j ) {
            residual += std::abs( u[cellIndexPE[j]]( 0, 0 ) - uOld[cellIndexPE[j]]( 0, 0 ) ) * _mesh->GetArea( cellIndexPE[j] );
        }
        MPI_Allreduce( &residual, &residualFull, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if( _settings->GetMyPE() == 0 ) {
            log->info( "{:03.8f}   {:01.5e}   {:01.5e}", t, residualFull, residualFull / dt );
            if( _settings->HasReferenceFile() && timeIndex % _settings->GetWriteFrequency() == 1 ) this->WriteErrors( refinementLevel );
        }
    }

    // write final error
    if( _settings->HasReferenceFile() ) this->WriteErrors( refinementLevel );

    // MPI Broadcast final moment vectors to all PEs
    for( unsigned j = 0; j < _nCells; ++j ) {
        MPI_Bcast( u[j].GetPointer(), int( _nStates * _nTotal ), MPI_DOUBLE, PEforCell[j], MPI_COMM_WORLD );
    }

    if( _settings->GetMyPE() != 0 ) return;

    // save final moments on uNew
    uNew = u;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    log->info( "" );
    log->info( "Finished!" );
    log->info( "" );
    log->info( "Runtime: {0}s", std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 );

    Matrix meanAndVar;
    Matrix meanAndVarErrors = Matrix( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
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
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k, nTotal[refinementLevel[j]] ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i, j ) += tmp[i] * phiTildeWf( k, 0 );
            }
        }

        // variance
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k, nTotal[refinementLevel[j]] ) );
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

    if( _settings->HasReferenceFile() ) {
        Vector a( 2 );
        a[0] = -0.05;
        a[1] = -0.5;
        Vector b( 2 );
        b[0]         = 1.05;
        b[1]         = 0.5;
        auto l1Error = this->CalculateError( meanAndVar, 1, a, b );
        auto l2Error = this->CalculateError( meanAndVar, 2, a, b );
        log->info( "\nExpectation Value error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = 0; i < _nStates; ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
        log->info( "\nVariance error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = _nStates; i < 2 * _nStates; ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
        meanAndVarErrors = this->CalculateErrorField( meanAndVar, 2 );

        _mesh->Export( meanAndVarErrors, "_errors" );
    }

    if( _settings->GetProblemType() == P_SHALLOWWATER_2D )
        _mesh->ExportShallowWater( meanAndVar );
    else {
        _mesh->Export( meanAndVar, "" );
    }

    // loop over all cells and check refinement indicator
    if( useAdaptivity ) {    // master determines refinement level for now
        for( unsigned j = 0; j < _nCells; ++j ) {
            // for( unsigned j = 0; j < static_cast<unsigned>( cellIndexPE.size() ); ++j ) {
            refinementIndicatorPlot( 0, j ) = std::fabs( u[j]( 0, _nTotal - 1 ) ) + std::fabs( u[j]( 0, _nTotal - 2 ) );    // modify for multiD
            refinementIndicatorPlot( 1, j ) = double( refinementLevel[j] );
        }
        _mesh->Export( refinementIndicatorPlot, "_refinementIndicator" );
    }

    this->Export( uNew, _lambda );

    unsigned evalCell = 2404;
    if( _settings->GetNumCells() > evalCell ) {
        std::ofstream outXi( "../results/xiGrid" );

        Vector xiEta( _settings->GetNDimXi() );
        unsigned n;

        unsigned plotState  = 0;
        unsigned nQFine     = 100;
        unsigned nQOriginal = _settings->GetNQuadPoints();
        _settings->SetNQuadPoints( nQFine );
        Closure* closurePlot = Closure::Create( _settings );
        Matrix testLambda( _nStates, _nTotal, 0.0 );
        testLambda( plotState, 1 ) = 1.0;
        _mesh->PlotInXi( closurePlot->U( closurePlot->EvaluateLambda( _lambda[evalCell] ) ), plotState );
        std::vector<Polynomial*> quad = closurePlot->GetQuadrature();
        testLambda                    = _lambda[evalCell];
        auto uPlot                    = closurePlot->U( closurePlot->EvaluateLambda( testLambda ) );

        for( unsigned k = 0; k < uPlot.columns(); ++k ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index = unsigned( ( k - k % unsigned( std::pow( nQFine, l ) ) ) / unsigned( std::pow( nQFine, l ) ) ) % nQFine;
                xiEta[l]       = quad[n]->GetNodes()[index];
                outXi << xiEta[l] << " ";
            }
            outXi << uPlot( 0, k ) << std::endl;
        }

        outXi.close();

        _settings->SetNQuadPoints( nQOriginal );
    }
}

MatVec MomentSolver::DetermineMoments( unsigned nTotal ) const {
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    if( _settings->HasRestartFile() ) {
        u = this->ImportPrevMoments( nTotal );

        // if cells are Dirichlet cells and we increase nTotal on restart, we should use moments from initial condition:
        // if we use old moments we have inexact BCs for new truncation order
        if( nTotal != _nTotal ) {
            MatVec uIC = this->SetupIC();
            for( unsigned j = 0; j < _nCells; ++j ) {
                if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) u[j] = uIC[j];
            }
        }
    }
    else {
        u = this->SetupIC();
    }
    return u;
}

Settings* MomentSolver::DeterminePreviousSettings() const {
    Settings* prevSettings;
    if( _settings->HasRestartFile() ) {
        prevSettings = this->ImportPrevSettings();
        prevSettings->SetNStates( _settings->GetNStates() );
        prevSettings->SetGamma( _settings->GetGamma() );
        prevSettings->SetClosureType( _settings->GetClosureType() );
    }
    else {
        prevSettings = _settings;
    }
    return prevSettings;
}

Closure* MomentSolver::DeterminePreviousClosure( Settings* prevSettings ) const {
    Closure* prevClosure;
    if( prevSettings->GetNMoments() != _settings->GetNMoments() || prevSettings->GetNQTotal() != _settings->GetNQTotal() ) {
        prevClosure = Closure::Create( prevSettings );
    }
    else {
        prevClosure = _closure;
    }
    return prevClosure;
}

void MomentSolver::SetDuals( Settings* prevSettings, Closure* prevClosure, MatVec& u ) {
    unsigned maxIterations = _closure->GetMaxIterations();

    if( _settings->LoadLambda() ) {
        _lambda = this->ImportPrevDuals( prevSettings->GetNTotal() );
    }
    else {
        // compute dual states for given moment vector
        _lambda = MatVec( _nCells, Matrix( _nStates, prevSettings->GetNTotal(), 0.0 ) );

        // compute first initial guess
        Vector ds( _nStates );
        Vector u0( _nStates );
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                u0[l] = u[j]( l, 0 );
            }
            _closure->DS( ds, u0 );
            for( unsigned l = 0; l < _nStates; ++l ) {
                _lambda[j]( l, 0 ) = ds[l];
            }
        }

        // Converge initial condition entropy variables for One Shot IPM or if truncation order is increased
        if( _settings->GetMaxIterations() == 1 || prevSettings->GetNMoments() != _settings->GetNMoments() ) {
            prevClosure->SetMaxIterations( 10000 );
            for( unsigned j = 0; j < _nCells; ++j ) {
                prevClosure->SolveClosureSafe( _lambda[j], u[j], prevSettings->GetNTotal(), prevSettings->GetNQTotal() );
            }
            prevClosure->SetMaxIterations( maxIterations );
        }
    }
    // for restart with increased number of moments reconstruct solution at finer quad points and compute moments for new truncation order
    if( prevSettings->GetNMoments() != _settings->GetNMoments() ) {
        MatVec uQFullProc = MatVec( _nCells, Matrix( _nStates, _settings->GetNQTotal() ) );
        if( maxIterations == 1 ) _closure->SetMaxIterations( 10000 );    // if one shot IPM is used, make sure that initial duals are converged
        prevSettings->SetNQuadPoints( _settings->GetNQuadPoints() );
        Closure* intermediateClosure = Closure::Create( prevSettings );    // closure with old nMoments and new Quadrature set
        for( unsigned j = 0; j < _nCells; ++j ) {
            _closure->U( uQFullProc[j], intermediateClosure->EvaluateLambda( _lambda[j] ) );    // solution at fine Quadrature nodes
            auto uCurrent = uQFullProc[j] * _closure->GetPhiTildeWf();
            u[j].resize( _nStates, _nTotal );
            u[j] = uCurrent;    // new Moments of size new nMoments

            // compute lambda with size newMoments
            Matrix lambdaOld = _lambda[j];
            _lambda[j].resize( _nStates, _nTotal );
            for( unsigned s = 0; s < _nStates; ++s ) {
                for( unsigned i = prevSettings->GetNTotal(); i < _settings->GetNTotal(); ++i ) {
                    _lambda[j]( s, i ) = 0.0;
                }
            }
            _closure->SolveClosureSafe( _lambda[j], u[j], _settings->GetNTotal(), _settings->GetNQTotal() );
        }
        _closure->SetMaxIterations( maxIterations );
        // delete reload closures and settings
        delete intermediateClosure;
        delete prevSettings;
    }
    if( prevSettings->GetNMoments() != _settings->GetNMoments() || prevSettings->GetNQTotal() != _settings->GetNQTotal() ) delete prevClosure;
}

void MomentSolver::numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n, unsigned nTotal ) {
    // out += _problem->G( u1, u2, nUnit, n ) * _closure->GetPhiTildeWf();
    multOnPENoReset( _problem->G( u1, u2, nUnit, n ), _closure->GetPhiTildeWf(), out, _settings->GetKStart(), _settings->GetKEnd(), nTotal );
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

MatVec MomentSolver::SetupIC() const {
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

void MomentSolver::Export( const MatVec& u, const MatVec& lambda ) const {
    std::shared_ptr<spdlog::logger> moment_writer = spdlog::get( "moments" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < _nTotal; ++k ) {
                line << std::setprecision( std::numeric_limits<double>::digits10 ) << u[i]( j, k ) << ",";
            }
        }
        moment_writer->info( line.str() );
    }
    moment_writer->flush();
    std::shared_ptr<spdlog::logger> dual_writer = spdlog::get( "duals" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < _nTotal; ++k ) {
                line << std::setprecision( std::numeric_limits<double>::digits10 ) << lambda[i]( j, k ) << ",";
            }
        }
        dual_writer->info( line.str() );
    }
    dual_writer->flush();
}

Settings* MomentSolver::ImportPrevSettings() const {
    auto file = std::ifstream( _settings->GetRestartFile() );
    std::string line;
    std::stringstream prevSettingsStream;
    bool configSection = false;
    while( std::getline( file, line ) ) {
        if( line.find( "Config file" ) != std::string::npos ) {
            configSection = true;
            for( unsigned i = 0; i < 4; ++i ) std::getline( file, line );
        }
        else if( configSection && line.find( "==================================" ) != std::string::npos ) {
            break;
        }
        if( configSection ) {
            line.erase( 0, line.find_first_of( '|' ) + 1 );
            prevSettingsStream << line << std::endl;
        }
    }
    std::istringstream inputStream( prevSettingsStream.str() );
    Settings* prevSettings = new Settings( inputStream );

    return prevSettings;
}

MatVec MomentSolver::ImportPrevMoments( unsigned nPrevTotal ) const {
    MatVec u( _nCells, Matrix( _nStates, nPrevTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_moments" );
    std::string line;
    for( unsigned j = 0; j < _nCells; ++j ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned i = 0; i < nPrevTotal; ++i ) {
                std::getline( lineStream, cell, ',' );
                u[j]( s, i ) = std::stod( cell );
            }
        }
    }
    file.close();
    return u;
}

MatVec MomentSolver::ImportPrevDuals( unsigned nPrevTotal ) {
    MatVec lambda( _nCells, Matrix( _nStates, nPrevTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_duals" );
    std::string line;
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < nPrevTotal; ++k ) {
                std::getline( lineStream, cell, ',' );
                lambda[i]( j, k ) = std::stod( cell );
            }
        }
    }
    file.close();
    return lambda;
}

Matrix MomentSolver::CalculateErrorField( const Matrix& solution, unsigned LNorm ) const {
    Matrix error( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    Vector refNorm( 2 * _nStates, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        switch( LNorm ) {
            case 1:
                for( unsigned s = 0; s < _nStates; ++s ) {
                    error( s, j )            = std::fabs( ( solution( s, j ) - _referenceSolution[j][s] ) );
                    error( _nStates + s, j ) = std::fabs( ( solution( _nStates + s, j ) - _referenceSolution[j][s + _nStates] ) );
                    refNorm[s] += std::fabs( _referenceSolution[j][s] ) * _mesh->GetArea( j );
                    refNorm[s + _nStates] += std::fabs( _referenceSolution[j][s + _nStates] ) * _mesh->GetArea( j );
                }
                break;
            case 2:
                for( unsigned s = 0; s < _nStates; ++s ) {
                    error( s, j )            = std::pow( ( solution( s, j ) - _referenceSolution[j][s] ), 2 );
                    error( _nStates + s, j ) = std::pow( ( solution( _nStates + s, j ) - _referenceSolution[j][s + _nStates] ), 2 );
                    refNorm[s] += std::pow( _referenceSolution[j][s], 2 ) * _mesh->GetArea( j );
                    refNorm[s + _nStates] += std::pow( _referenceSolution[j][s + _nStates], 2 ) * _mesh->GetArea( j );
                }
                break;
            default: exit( EXIT_FAILURE );
        }
    }

    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned s = 0; s < _nStates; ++s ) {
            error( s, j )            = error( s, j ) / refNorm[s];
            error( _nStates + s, j ) = error( _nStates + s, j ) / refNorm[s + _nStates];
        }
    }

    return error;
}

Vector MomentSolver::CalculateError( const Matrix& solution, unsigned LNorm, const Vector& a, const Vector& b ) const {
    Vector error( 2 * _nStates, 0.0 );
    Vector refNorm( 2 * _nStates, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _mesh->GetGrid()[j]->GetCenter()[0] > a[0] && _mesh->GetGrid()[j]->GetCenter()[0] < b[0] && _mesh->GetGrid()[j]->GetCenter()[1] > a[1] &&
            _mesh->GetGrid()[j]->GetCenter()[1] < b[1] ) {

            switch( LNorm ) {
                case 1:
                    for( unsigned s = 0; s < 2 * _nStates; ++s ) {
                        error[s] += std::fabs( ( solution( s, j ) - _referenceSolution[j][s] ) ) * _mesh->GetArea( j );
                        refNorm[s] += std::fabs( _referenceSolution[j][s] ) * _mesh->GetArea( j );
                    }
                    break;
                case 2:
                    for( unsigned s = 0; s < 2 * _nStates; ++s ) {
                        error[s] += std::pow( ( solution( s, j ) - _referenceSolution[j][s] ), 2 ) * _mesh->GetArea( j );
                        refNorm[s] += std::pow( _referenceSolution[j][s], 2 ) * _mesh->GetArea( j );
                    }

                    break;
                default: exit( EXIT_FAILURE );
            }
        }
    }
    for( unsigned s = 0; s < 2 * _nStates; ++s ) {
        error[s] = std::pow( error[s] / refNorm[s], 1.0 / double( LNorm ) );
    }
    return error;
}

void MomentSolver::WriteErrors( const VectorU& refinementLevel ) {
    // define rectangle for error computation
    Vector a( 2 );
    a[0] = -0.05;
    a[1] = -0.5;
    Vector b( 2 );
    b[0]                = 1.05;
    b[1]                = 0.5;
    auto l1ErrorMeanLog = spdlog::get( "l1ErrorMean" );
    auto l2ErrorMeanLog = spdlog::get( "l2ErrorMean" );
    // auto lInfErrorMeanLog = spdlog::get( "lInfErrorMean" );
    auto l1ErrorVarLog = spdlog::get( "l1ErrorVar" );
    auto l2ErrorVarLog = spdlog::get( "l2ErrorVar" );
    // auto lInfErrorVarLog  = spdlog::get( "lInfErrorVar" );

    Matrix meanAndVar = Matrix( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    Matrix phiTildeWf = _closure->GetPhiTildeWf();
    Vector tmp( _nStates, 0.0 );
    VectorU nTotal = _settings->GetNTotalRefinementLevel();
    for( unsigned j = 0; j < _nCells; ++j ) {
        // expected value
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k, nTotal[refinementLevel[j]] ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i, j ) += tmp[i] * phiTildeWf( k, 0 );
            }
        }

        // variance
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k, nTotal[refinementLevel[j]] ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i + _nStates, j ) += pow( tmp[i] - meanAndVar( i, j ), 2 ) * phiTildeWf( k, 0 );
            }
        }
    }
    auto l1Error = this->CalculateError( meanAndVar, 1, a, b );
    auto l2Error = this->CalculateError( meanAndVar, 2, a, b );
    // auto lInfError = this->CalculateError( meanAndVar, 0, a, b );

    std::ostringstream osL1ErrorMean, osL2ErrorMean, osLInfErrorMean, osL1ErrorVar, osL2ErrorVar, osLInfErrorVar;
    for( unsigned i = 0; i < _nStates; ++i ) {
        osL1ErrorMean << std::scientific << l1Error[i] << "\t";
        osL2ErrorMean << std::scientific << l2Error[i] << "\t";
        // osLInfErrorMean << std::scientific << lInfError[i] << "\t";
        osL1ErrorVar << std::scientific << l1Error[i + _nStates] << "\t";
        osL2ErrorVar << std::scientific << l2Error[i + _nStates] << "\t";
        // osLInfErrorVar << std::scientific << lInfError[i + _nStates] << "\t";
    }

    l1ErrorMeanLog->info( osL1ErrorMean.str() );
    l2ErrorMeanLog->info( osL2ErrorMean.str() );
    // lInfErrorMeanLog->info( osLInfErrorMean.str() );
    l1ErrorVarLog->info( osL1ErrorVar.str() );
    l2ErrorVarLog->info( osL2ErrorVar.str() );
    // lInfErrorVarLog->info( osLInfErrorVar.str() );
}
