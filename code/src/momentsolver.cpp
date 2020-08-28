#include "momentsolver.h"
#include <mpi.h>

MomentSolver::MomentSolver( Settings* settings, Mesh* mesh, Problem* problem ) : _settings( settings ), _mesh( mesh ), _problem( problem ) {
    _log         = spdlog::get( "event" );
    _nCells      = _settings->GetNumCells();
    _tStart      = 0.0;
    _tEnd        = _settings->GetTEnd();
    _nStates     = _settings->GetNStates();
    _nQuadPoints = _settings->GetNQuadPoints();

    _nTotal         = _settings->GetNTotal();
    _nTotalForRef   = _settings->GetNTotalRefinementLevel();
    _nQTotalForRef  = _settings->GetNQTotalForRef();
    _nMultiElements = _settings->GetNMultiElements();

    _cellIndexPE = _settings->GetCellIndexPE();

    _closure = Closure::Create( _settings );
    _nQTotal = _settings->GetNQTotal();
    _time    = TimeSolver::Create( _settings, _mesh, _problem );

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
    unsigned retCounter = 0;    // counter for retardation level
    // unsigned retLevel        = _settings->GetResidualRetardation( 0 );
    VectorU refinementLevel( _nCells, _settings->GetNRefinementLevels( retCounter ) - 1 );       // vector carries refinement level for each cell
    VectorU refinementLevelOld( _nCells, _settings->GetNRefinementLevels( retCounter ) - 1 );    // vector carries old refinement level for each cell
    VectorU refinementLevelTransition( _nCells, _settings->GetNRefinementLevels( retCounter ) - 1 );

    // set Dirichlet Cells to finest refinement level
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
            refinementLevel[j]           = _settings->GetNRefinementLevels() - 1;
            refinementLevelOld[j]        = _settings->GetNRefinementLevels() - 1;
            refinementLevelTransition[j] = _settings->GetNRefinementLevels() - 1;
        }
    }

    auto log                                  = spdlog::get( "event" );
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    // create settings+closure for restart as well as solution fields
    Settings* prevSettings = DeterminePreviousSettings();
    Closure* prevClosure   = DeterminePreviousClosure( prevSettings );
    if( _settings->HasRestartFile() ) _tStart = prevSettings->GetTEnd();
    MatTens u = DetermineMoments( prevSettings->GetNTotal() );
    SetDuals( prevSettings, prevClosure, u );
    MatTens uQ    = MatTens( _nCells, Tensor( _nStates, _nMultiElements, _settings->GetNqPE() ) );
    MatTens uQNew = MatTens( _nCells, Tensor( _nStates, _nMultiElements, _settings->GetNqPE() ) );
    // slopes for slope limiter
    MatTens duQx( _nCells, Tensor( _nStates, _nMultiElements, _settings->GetNqPE(), 0.0 ) );
    MatTens duQy( _nCells, Tensor( _nStates, _nMultiElements, _settings->GetNqPE(), 0.0 ) );

    std::vector<int> PEforCell = _settings->GetPEforCell();

    MatTens uNew = u;
    MatTens uOld = u;

    // set up function pointer for right hand side
    auto numFluxPtr = std::bind( &MomentSolver::numFlux, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

    if( _settings->GetMyPE() == 0 ) log->info( "{:10}   {:10}", "t", "residual" );

    // init time and residual
    double t           = _tStart;
    unsigned timeIndex = 0;
    double dt;
    double minResidual  = _settings->GetMinResidual();
    double residualFull = minResidual + 1.0;

    // perform initial step for regularization
    if( _settings->HasRegularization() ) PerformInitialStep( refinementLevel, uNew );

    // Begin time loop
    while( t < _tEnd && residualFull > minResidual ) {

        double residual = 0;
        // Solve dual problem
#pragma omp parallel for schedule( dynamic, 10 )
        for( unsigned j = 0; j < static_cast<unsigned>( _cellIndexPE.size() ); ++j ) {
            if( _mesh->GetBoundaryType( _cellIndexPE[j] ) == BoundaryType::DIRICHLET && timeIndex > 0 && !_settings->HasSource() ) continue;
            // std::cout << "Cell " << _cellIndexPE[j] << ", lambda = " << _lambda[_cellIndexPE[j]] << ", u = " << u[_cellIndexPE[j]] << std::endl;
            _closure->SolveClosureSafe( _lambda[_cellIndexPE[j]], u[_cellIndexPE[j]], refinementLevel[_cellIndexPE[j]] );
            // std::cout << "result = " << _lambda[_cellIndexPE[j]] << std::endl;
        }

        // MPI Broadcast lambdas to all PEs
        for( unsigned j = 0; j < _nCells; ++j ) {
            uOld[j] = u[j];    // save old Moments for residual computation
            MPI_Bcast( _lambda[j].GetPointer(), int( _nStates * _nTotal * _nMultiElements ), MPI_DOUBLE, PEforCell[j], MPI_COMM_WORLD );
        }

        // save old refinement levels
        refinementLevelOld = refinementLevel;

        // determine refinement level of cells on current PE
        for( unsigned j = 0; j < static_cast<unsigned>( _cellIndexPE.size() ); ++j ) {
            // if( _mesh->GetBoundaryType( _cellIndexPE[j] ) == BoundaryType::DIRICHLET && timeIndex > 0 ) continue;
            double indicator;
            if( _settings->GetNTotal() == 1 )    // turn off adaptivity if just one moment is used
                indicator = 0.0;
            else
                indicator = ComputeRefIndicator( u[_cellIndexPE[j]], refinementLevel[_cellIndexPE[j]] );

            if( indicator > _settings->GetRefinementThreshold() &&
                refinementLevel[_cellIndexPE[j]] < _settings->GetNRefinementLevels( retCounter ) - 1 &&
                _mesh->GetBoundaryType( _cellIndexPE[j] ) != BoundaryType::DIRICHLET )
                refinementLevel[_cellIndexPE[j]] += 1;
            else if( indicator < _settings->GetCoarsenThreshold() && refinementLevel[_cellIndexPE[j]] > 0 &&
                     _mesh->GetBoundaryType( _cellIndexPE[j] ) != BoundaryType::DIRICHLET )
                refinementLevel[_cellIndexPE[j]] -= 1;
        }

        // broadcast refinemt level to all PEs
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Bcast( &refinementLevel[j], 1, MPI_UNSIGNED, PEforCell[j], MPI_COMM_WORLD );
        }

        // determine transition refinement level to ensure that neighboring cells of high refinement level have fine uQ reconstruction
        // important for FV stencil in Advance function
        for( unsigned j = 0; j < static_cast<unsigned>( _cellIndexPE.size() ); ++j ) {
            refinementLevelTransition[_cellIndexPE[j]] = refinementLevel[_cellIndexPE[j]];
            auto neighbors                             = _mesh->GetNeighborIDs( _cellIndexPE[j] );
            unsigned maxRefLevelNghs                   = refinementLevel[_cellIndexPE[j]];
            for( unsigned ngh = 0; ngh < neighbors.size(); ++ngh ) {
                if( neighbors[ngh] != _nCells ) {
                    if( maxRefLevelNghs < refinementLevel[neighbors[ngh]] ) {
                        maxRefLevelNghs                            = refinementLevel[neighbors[ngh]];
                        refinementLevelTransition[_cellIndexPE[j]] = maxRefLevelNghs;
                    }
                }
            }
        }
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Bcast( &refinementLevelTransition[j], 1, MPI_UNSIGNED, PEforCell[j], MPI_COMM_WORLD );
        }

        // compute solution at quad points
        for( unsigned j = 0; j < _nCells; ++j ) {
            uQ[j] = _closure->U( _closure->EvaluateLambdaOnPE( _lambda[j], refinementLevelOld[j], refinementLevelTransition[j] ) );
        }

        // compute partial moment vectors on each PE (for inexact dual variables)
        for( unsigned j = 0; j < _nCells; ++j ) {
            // recompute moments with inexact dual variables
            u[j] = uQ[j] * _closure->GetPhiTildeWfAtRef( refinementLevel[j] );
        }

        // determine time step size
        double dtMinOnPE = 1e10;
        double dtCurrent;
        for( unsigned j = 0; j < _nCells; ++j ) {
            dtCurrent = _problem->ComputeDt( uQ[j], _mesh->GetMaxEdge( j ) / _mesh->GetArea( j ), refinementLevel[j] );
            if( dtCurrent < dtMinOnPE ) dtMinOnPE = dtCurrent;
        }
        MPI_Allreduce( &dtMinOnPE, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
        _settings->SetDT( dt );
        t += dt;
        ++timeIndex;

        // perform time update flux
        _time->Advance( numFluxPtr, uQNew, u, uQ, dt, refinementLevel );

        // perform time update source
        if( _settings->HasSource() ) {
            this->Source( uQNew, uQ, dt, t, refinementLevel );
        }

        // compute partial moments from time updated solution at quad points on each PE
        for( unsigned j = 0; j < _nCells; ++j ) {
            uNew[j] = uQNew[j] * _closure->GetPhiTildeWfAtRef( refinementLevel[j] );
        }

        // perform reduction onto u to obtain full moments on PE PEforCell[j], which is the PE that solves the dual problem for cell j
        for( unsigned j = 0; j < _nCells; ++j ) {
            MPI_Reduce( uNew[j].GetPointer(),
                        u[j].GetPointer(),
                        int( _nStates * _nMultiElements * _nTotal ),
                        MPI_DOUBLE,
                        MPI_SUM,
                        PEforCell[j],
                        MPI_COMM_WORLD );
        }

        if( _settings->HasRegularization() ) {
            // add eta*lambda to obtain old moments
            for( unsigned j = 0; j < static_cast<unsigned>( _cellIndexPE.size() ); ++j ) {
                u[_cellIndexPE[j]] = u[_cellIndexPE[j]].Add( _lambda[_cellIndexPE[j]] * _settings->GetRegularizationStrength(),
                                                             _nStates,
                                                             _nMultiElements,
                                                             _nTotalForRef[refinementLevelOld[_cellIndexPE[j]]] );
            }
        }

        // compute residual
        for( unsigned j = 0; j < _cellIndexPE.size(); ++j ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                residual += std::abs( u[_cellIndexPE[j]]( 0, l, 0 ) - uOld[_cellIndexPE[j]]( 0, l, 0 ) ) * _mesh->GetArea( _cellIndexPE[j] );
            }
        }

        MPI_Allreduce( &residual, &residualFull, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        if( _settings->GetMyPE() == 0 ) {
            log->info( "{:03.8f}   {:01.5e}   {:01.5e}", t, residualFull, residualFull / dt );
            if( _settings->HasReferenceFile() && timeIndex % _settings->GetWriteFrequency() == 1 ) this->WriteErrors( refinementLevel, t );
            if( false && _settings->WriteInTime() && timeIndex % _settings->GetWriteFrequency() == 1 ) {
                Matrix meanAndVar = WriteMeanAndVar( refinementLevel, t, true, true );
                _mesh->Export( meanAndVar, "_" + std::to_string( timeIndex ) );
                if( _settings->GetNRefinementLevels() > 1 ) ExportRefinementIndicator( refinementLevel, u, timeIndex );
            }
        }

        // if current retardation level fulfills residual condition, then increase retardation level
        if( residualFull < _settings->GetResidualRetardation( retCounter ) && retCounter < _settings->GetNRetardationLevels() - 1 ) retCounter += 1;
    }

    // write final error
    if( _settings->HasReferenceFile() ) this->WriteErrors( refinementLevel, t );

    // MPI Broadcast final moment vectors to all PEs
    for( unsigned j = 0; j < _nCells; ++j ) {
        MPI_Bcast( u[j].GetPointer(), int( _nStates * _nMultiElements * _nTotal ), MPI_DOUBLE, PEforCell[j], MPI_COMM_WORLD );
    }

    if( _settings->GetMyPE() != 0 ) return;

    // save final moments on uNew
    uNew = u;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    log->info( "" );
    log->info( "Finished!" );
    log->info( "" );
    log->info( "Runtime: {0}s", std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 );

    // compute mean and variance numerical + exact (if exact solution specified)
    Matrix meanAndVar = WriteMeanAndVar( refinementLevel, t, true, false );
    // write exact solution on reference field
    if( _settings->HasExactSolution() && !_settings->HasReferenceFile() ) {
        _referenceSolution.resize( _nCells );
        for( unsigned j = 0; j < _nCells; ++j ) {
            _referenceSolution[j].resize( 2 * _nStates );
            for( unsigned s = 2 * _nStates; s < 4 * _nStates; ++s ) {
                _referenceSolution[j][s - 2 * _nStates] = meanAndVar( s, j );
            }
        }
    }

    if( _settings->HasReferenceFile() || _settings->HasExactSolution() ) {
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
        Matrix meanAndVarErrors = this->CalculateErrorField( meanAndVar, 2 );

        // export error plot
        _mesh->Export( meanAndVarErrors, "_errors" );

        // export 2nd derivative error
        if( _settings->HasReferenceFile() ) Write2ndDerMeanAndVar( refinementLevel, t );
    }

    // export mean and variance
    if( _settings->GetProblemType() == P_SHALLOWWATER_2D )
        _mesh->ExportShallowWater( meanAndVar );
    else {
        Matrix meanVarT = WriteMeanAndVar( refinementLevel, t, true, true );
        _mesh->Export( meanVarT, "" );
    }

    // WriteGradientsScalarField( meanAndVar );

    // export solution fields
    this->Export( uNew, _lambda, refinementLevel );

    // export refinement indicator
    if( _settings->GetNRefinementLevels() > 1 ) {
        ExportRefinementIndicator( refinementLevel, u, 1000 );
    }

    unsigned evalCell = 333;
    if( _settings->GetNumCells() > evalCell ) {
        std::ofstream outXi( "../results/xiGrid" );

        Vector xiEta( _settings->GetNDimXi() );
        unsigned n;

        unsigned plotState  = 0;
        unsigned nQFine     = 1000;
        unsigned nQOriginal = _settings->GetNQuadPoints();
        _settings->SetNQuadPoints( nQFine );
        Closure* closurePlot = Closure::Create( _settings );
        Tensor testLambda( _nStates, _nMultiElements, _nTotal, 0.0 );
        for( unsigned l = 0; l < _nMultiElements; ++l ) testLambda( plotState, l, 1 ) = 1.0;
        _mesh->PlotInXi( closurePlot->U( closurePlot->EvaluateLambda( _lambda[evalCell] ) ), plotState );
        std::vector<Polynomial*> quad = closurePlot->GetQuadrature();
        testLambda                    = _lambda[evalCell];
        auto uPlot                    = closurePlot->U( closurePlot->EvaluateLambda( testLambda ) );

        for( unsigned k = 0; k < uPlot.columns(); ++k ) {
            for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index = unsigned( ( k - k % unsigned( std::pow( nQFine, m ) ) ) / unsigned( std::pow( nQFine, m ) ) ) % nQFine;
                xiEta[m]       = quad[n]->GetNodes()[index];
                outXi << xiEta[m] << " ";
            }
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                outXi << uPlot( 0, l, k ) << std::endl;
            }
        }

        outXi.close();
        delete closurePlot;

        _settings->SetNQuadPoints( nQOriginal );
    }
}

void MomentSolver::Source( MatTens& uQNew, const MatTens& uQ, double dt, double t, const VectorU& refLevel ) const {
    unused( uQ );

    //#pragma omp parallel for
    auto uQTilde = uQNew;

    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) continue;
        Tensor out = _problem->Source( uQNew[j], _mesh->GetCenterPos( j ), t, refLevel[j] );    //  use uQ or uQNew?
        uQNew[j]   = uQNew[j] + out * dt;
    }
}

void MomentSolver::numFlux( Matrix& out, const Matrix& g, unsigned level ) {
    unused( level );

    out += g;
}

double MomentSolver::ComputeRefIndicator( const Tensor& u, unsigned refLevel ) const {
    double indicator = 0;
    double P         = 1.0 / _nMultiElements;
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        if( _settings->GetNDimXi() == 1 ) {
            if( _nTotalForRef[refLevel] - 2 > 0 ) {    // take last two moments as indicator to account for even/odd effects
                indicator += std::fabs( u( 0, l, _nTotalForRef[refLevel] - 1 ) ) + std::fabs( u( 0, l, _nTotalForRef[refLevel] - 2 ) );
            }
            else if( _nTotalForRef[refLevel] - 1 > 0 ) {    // if only moment up to order 1 exists, take this twice
                indicator += 2.0 * std::fabs( u( 0, l, _nTotalForRef[refLevel] - 1 ) );
            }
        }
        else {    // in multi-D, take moments from highest level
            unsigned prevRefinementLevel;
            if( refLevel == 0 ) {
                prevRefinementLevel = 1;
            }
            else {
                prevRefinementLevel = _nTotalForRef[refLevel - 1];
            }
            for( unsigned i = prevRefinementLevel; i < _nTotalForRef[refLevel]; ++i ) {
                indicator += std::fabs( u( 0, l, i ) );
            }
        }
    }
    return indicator * P;
}

void MomentSolver::ExportRefinementIndicator( const VectorU& refinementLevel, const MatTens& u, unsigned index ) const {
    unused( u );

    // loop over all cells and check refinement indicator
    Matrix refinementIndicatorPlot( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        double indicator                = 1.0;
        refinementIndicatorPlot( 0, j ) = indicator;    // modify for multiD
        refinementIndicatorPlot( 1, j ) = double( refinementLevel[j] );
    }
    _mesh->Export( refinementIndicatorPlot, "_refinementIndicator" + std::to_string( index ) );
}

Matrix MomentSolver::WriteMeanAndVar( const VectorU& refinementLevel, double t, bool writeExact, bool writeTemperature ) const {
    Matrix meanAndVar;

    // write out temperature if Euler2D
    unsigned nStates;
    if( _settings->GetProblemType() == P_EULER_2D && writeTemperature )
        nStates = _nStates + 1;
    else
        nStates = _nStates;
    std::vector<std::vector<unsigned>> indicesQ;
    if( _settings->HasExactSolution() && writeExact ) {
        meanAndVar = Matrix( 4 * nStates, _mesh->GetNumCells(), 0.0 );
    }
    else {
        meanAndVar = Matrix( 2 * nStates, _mesh->GetNumCells(), 0.0 );
    }
    Vector tmp( _nStates, 0.0 );
    for( unsigned j = 0; j < _nCells; ++j ) {
        Matrix phiTildeWf = _closure->GetPhiTildeWfAtRef( refinementLevel[j], true );
        double P          = 1.0 / _nMultiElements;
        // expected value
        for( unsigned k = 0; k < _settings->GetNQTotalForRef( refinementLevel[j] ); ++k ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], l, k, _nTotalForRef[refinementLevel[j]] ) );
                for( unsigned i = 0; i < _nStates; ++i ) {
                    meanAndVar( i, j ) += tmp[i] * phiTildeWf( k, 0 ) * P;
                }
                if( _settings->GetProblemType() == P_EULER_2D && writeTemperature ) {
                    double rho  = tmp[0];
                    double rhoU = tmp[1];
                    double rhoV = tmp[2];
                    double rhoE = tmp[3];
                    double p    = 0.4 * ( rhoE - 0.5 * ( pow( rhoU, 2 ) + pow( rhoV, 2 ) ) / rho );
                    double T    = p / ( rho * 287.87 );
                    meanAndVar( nStates - 1, j ) += T * phiTildeWf( k, 0 ) * P;
                }
            }
        }

        // variance
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned k = 0; k < _settings->GetNQTotalForRef( refinementLevel[j] ); ++k ) {
                _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], l, k, _nTotalForRef[refinementLevel[j]] ) );
                for( unsigned i = 0; i < _nStates; ++i ) {
                    meanAndVar( i + nStates, j ) += pow( tmp[i] - meanAndVar( i, j ), 2 ) * phiTildeWf( k, 0 ) * P;
                }
                if( _settings->GetProblemType() == P_EULER_2D && writeTemperature ) {
                    double rho  = tmp[0];
                    double rhoU = tmp[1];
                    double rhoV = tmp[2];
                    double rhoE = tmp[3];
                    double p    = 0.4 * ( rhoE - 0.5 * ( pow( rhoU, 2 ) + pow( rhoV, 2 ) ) / rho );
                    double T    = p / ( rho * 287.87 );
                    meanAndVar( 2 * nStates - 1, j ) += pow( T - meanAndVar( nStates - 1, j ), 2 ) * phiTildeWf( k, 0 ) * P;
                }
            }
        }
    }
    if( _settings->HasExactSolution() && writeExact ) {
        // store xGrid on vector
        Matrix xGrid( _nCells, _settings->GetMeshDimension() );
        for( unsigned j = 0; j < _nCells; ++j ) {
            Vector midPj = _mesh->GetCenterPos( j );
            for( unsigned s = 0; s < _settings->GetMeshDimension(); ++s ) {
                xGrid( j, s ) = midPj[s];
            }
        }
        Vector xiEta( _settings->GetNDimXi() );
        std::vector<Polynomial*> quad( 2 );
        unsigned nQuadFine = 5 * _settings->GetNQuadPoints();    // define fine quadrature for exact solution
        quad[0]            = Polynomial::Create( _settings, nQuadFine, DistributionType::D_LEGENDRE );
        quad[1]            = Polynomial::Create( _settings, nQuadFine, DistributionType::D_HERMITE );
        unsigned n;

        // compute indices for quad points
        unsigned nQTotal = pow( nQuadFine, _settings->GetNDimXi() );
        indicesQ.resize( nQTotal );
        for( unsigned k = 0; k < nQTotal; ++k ) {
            indicesQ[k].resize( _settings->GetNDimXi() );
            for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                indicesQ[k][m] = unsigned( ( k - k % unsigned( std::pow( nQuadFine, m ) ) ) / unsigned( std::pow( nQuadFine, m ) ) ) % nQuadFine;
            }
        }

        for( unsigned k = 0; k < nQTotal; ++k ) {
            for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                xiEta[m] = quad[n]->GetNodes()[indicesQ[k][m]];
            }
            Matrix exactSolOnMesh = _problem->ExactSolution( t, xGrid, xiEta );
            for( unsigned j = 0; j < _nCells; ++j ) {
                // expected value exact
                for( unsigned i = 0; i < _nStates; ++i ) {
                    double wfXi = 1.0;
                    for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                        if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                        if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                        wfXi *= quad[n]->fXi( quad[n]->GetNodes()[indicesQ[k][m]] ) * quad[n]->GetWeights()[indicesQ[k][m]];
                    }

                    meanAndVar( 2 * nStates + i, j ) += exactSolOnMesh( j, i ) * wfXi;
                }
            }
        }
        // variance exact
        for( unsigned k = 0; k < nQTotal; ++k ) {
            for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                xiEta[m] = quad[n]->GetNodes()[indicesQ[k][m]];
            }

            Matrix exactSolOnMesh = _problem->ExactSolution( t, xGrid, xiEta );
            for( unsigned j = 0; j < _nCells; ++j ) {
                // expected value exact
                for( unsigned i = 0; i < _nStates; ++i ) {
                    double wfXi = 1.0;
                    for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                        if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                        if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                        wfXi *= quad[n]->fXi( quad[n]->GetNodes()[indicesQ[k][m]] ) * quad[n]->GetWeights()[indicesQ[k][m]];
                    }
                    meanAndVar( 3 * nStates + i, j ) += pow( exactSolOnMesh( j, i ) - meanAndVar( 2 * nStates + i, j ), 2 ) * wfXi;
                }
            }
        }
        delete quad[0];
        delete quad[1];
    }
    return meanAndVar;
}

MatTens MomentSolver::DetermineMoments( unsigned nTotal ) const {
    MatTens u( _nCells, Tensor( _nStates, _nMultiElements, _nTotal ) );
    if( _settings->HasRestartFile() ) {
        u = this->ImportPrevMoments( nTotal );

        // if cells are Dirichlet cells and we increase nTotal on restart, we should use moments from initial condition:
        // if we use old moments we have inexact BCs for new truncation order
        if( nTotal != _nTotal ) {
            MatTens uIC = this->SetupIC();
            for( unsigned j = 0; j < _nCells; ++j ) {
                if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                    u[j].resize( _nStates, _nMultiElements, _nTotal );
                    u[j] = uIC[j];
                }
            }
        }
    }
    else {
        u = this->SetupIC();
    }
    return u;
}

// TODO: refLevel is stored for all x <-> MPI ?
void MomentSolver::DetermineGradients( MatTens& duQx, MatTens& duQy, const MatTens& uQ, const VectorU& refLevel ) const {
    for( unsigned j = 0; j < _nCells; ++j ) {
        duQx[j].reset();
        duQy[j].reset();
        auto cell        = _mesh->GetCell( j );
        auto neighborIDs = cell->GetNeighborIDs();
        // compute derivative for every quadrature point on PE
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned n = 0; n < _nMultiElements; ++n ) {
                for( unsigned k = 0; k < _nQTotalForRef[refLevel[j]]; ++k ) {
                    // use all neighboring cells for stencil
                    for( unsigned idx = 0; idx < neighborIDs.size(); ++idx ) {
                        duQx[j]( s, n, k ) = duQx[j]( s, n, k ) + 0.5 * ( uQ[j]( s, n, k ) + uQ[neighborIDs[idx]]( s, n, k ) ) *
                                                                      cell->GetNormal( idx )[0] / cell->GetArea();
                        duQy[j]( s, n, k ) = duQy[j]( s, n, k ) + 0.5 * ( uQ[j]( s, n, k ) + uQ[neighborIDs[idx]]( s, n, k ) ) *
                                                                      cell->GetNormal( idx )[1] / cell->GetArea();
                    }
                }
            }
        }
    }
}

void MomentSolver::DetermineGradientsScalarField( Matrix& dux, Matrix& duy, const Matrix& u ) const {
    dux.reset();
    duy.reset();
    unsigned nStates = u.rows();
    //    unsigned numCells = _mesh->GetNumCells();

    for( unsigned s = 0; s < nStates; ++s ) {
        for( unsigned j = 0; j < _nCells; ++j ) {
            auto cell = _mesh->GetCell( j );
            if( cell->IsBoundaryCell() ) continue;
            auto neighborIDs = cell->GetNeighborIDs();
            // use all neighboring cells for stencil
            for( unsigned idx = 0; idx < neighborIDs.size(); ++idx ) {
                dux( s, j ) = dux( s, j ) + 0.5 * ( u( s, j ) + u( s, neighborIDs[idx] ) ) * cell->GetNormal( idx )[0] / cell->GetArea();
                duy( s, j ) = duy( s, j ) + 0.5 * ( u( s, j ) + u( s, neighborIDs[idx] ) ) * cell->GetNormal( idx )[1] / cell->GetArea();
            }
        }
    }
}

void MomentSolver::WriteGradientsScalarField( const Matrix& u ) const {
    unsigned nStates = u.rows();
    Matrix dux( nStates, _nCells, 0.0 );
    Matrix duy( nStates, _nCells, 0.0 );
    DetermineGradientsScalarField( dux, duy, u );
    _mesh->Export( dux, "_xDer" );
    _mesh->Export( duy, "_yDer" );
}

void MomentSolver::Write2ndDerMeanAndVar( const VectorU& refinementLevel, double t ) const {
    Matrix MeanVar( 2 * _nStates, _nCells );
    Matrix MeanVarExact( 2 * _nStates, _nCells );
    Matrix meanAndVar = WriteMeanAndVar( refinementLevel, t, true, false );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned j = 0; j < _nCells; ++j ) {
            MeanVar( s, j )                 = meanAndVar( 0 * _nStates + s, j );
            MeanVar( _nStates + s, j )      = meanAndVar( 1 * _nStates + s, j );
            MeanVarExact( s, j )            = _referenceSolution[j][s];
            MeanVarExact( _nStates + s, j ) = _referenceSolution[j][s + _nStates];
        }
    }

    Matrix dux( 2 * _nStates, _nCells, 0.0 );
    Matrix duy( 2 * _nStates, _nCells, 0.0 );
    DetermineGradientsScalarField( dux, duy, MeanVarExact - MeanVar );
    //_mesh->Export( dux, "_ErrorXDer" );
    //_mesh->Export( duy, "_ErrorYDer" );
    Matrix duxx( 2 * _nStates, _nCells, 0.0 );
    Matrix duyy( 2 * _nStates, _nCells, 0.0 );
    DetermineGradientsScalarField( duxx, duyy, dux );
    //_mesh->Export( duxx, "_ErrorXXDer" );
    //_mesh->Export( duyy, "_ErrorXYDer" );
    DetermineGradientsScalarField( duxx, duyy, duy );
    //_mesh->Export( duxx, "_ErrorYXDer" );
    //_mesh->Export( duyy, "_ErrorYYDer" );

    Vector l1Error( 2 * _nStates, 0.0 );
    Vector l2Error( 2 * _nStates, 0.0 );
    for( unsigned s = 0; s < 2 * _nStates; ++s ) {
        for( unsigned j = 0; j < _nCells; ++j ) {
            l1Error[s] += ( std::fabs( duxx( s, j ) ) + std::fabs( duyy( s, j ) ) ) * _mesh->GetArea( j );
            l2Error[s] += ( std::pow( duxx( s, j ), 2 ) + std::pow( duyy( s, j ), 2 ) ) * _mesh->GetArea( j );
        }
        l2Error[s] = sqrt( l2Error[s] );
    }

    _log->info( "\nExpectation 2nd Der error w.r.t reference solution:" );
    _log->info( "State   L1-error      L2-error" );
    for( unsigned i = 0; i < _nStates; ++i ) {
        _log->info( "{:1d}       {:01.9e}   {:01.9e}", i, l1Error[i], l2Error[i] );
    }
    _log->info( "\nVariance 2nd Der error w.r.t reference solution:" );
    _log->info( "State   L1-error      L2-error" );
    for( unsigned i = _nStates; i < 2 * _nStates; ++i ) {
        _log->info( "{:1d}       {:01.9e}   {:01.9e}", i, l1Error[i], l2Error[i] );
    }
}

Settings* MomentSolver::DeterminePreviousSettings() const {
    Settings* prevSettings;
    if( _settings->HasRestartFile() ) {
        prevSettings = this->ImportPrevSettings();
        prevSettings->SetNStates( _settings->GetNStates() );
        prevSettings->SetGamma( _settings->GetGamma() );
        prevSettings->SetClosureType( _settings->GetClosureType() );
        prevSettings->SetNumCells( _settings->GetNumCells() );
        prevSettings->SetInputFile( _settings->GetInputFile() );
        prevSettings->SetDT( _dt );
    }
    else {
        prevSettings = _settings;
    }
    return prevSettings;
}

Closure* MomentSolver::DeterminePreviousClosure( Settings* prevSettings ) const {
    Closure* prevClosure;
    if( prevSettings->GetMaxDegree() != _settings->GetMaxDegree() || prevSettings->GetNQTotal() != _settings->GetNQTotal() ) {
        prevClosure = Closure::Create( prevSettings );
    }
    else {
        prevClosure = _closure;
    }
    return prevClosure;
}

void MomentSolver::SetDuals( Settings* prevSettings, Closure* prevClosure, MatTens& u ) {
    unsigned maxIterations = _closure->GetMaxIterations();

    // load or compute duals of read-in solution
    if( _settings->LoadLambda() ) {
        _lambda = this->ImportPrevDuals( prevSettings->GetNTotal() );
    }
    else {
        // compute dual states for given moment vector
        _lambda = MatTens( _nCells, Tensor( _nStates, _nMultiElements, prevSettings->GetNTotal(), 0.0 ) );

        // compute first initial guess
        Vector ds( _nStates );
        Vector u0( _nStates );
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned n = 0; n < _nMultiElements; ++n ) {
                for( unsigned s = 0; s < _nStates; ++s ) {
                    u0[s] = u[j]( s, n, 0 );
                }
                _closure->DS( ds, u0 );
                for( unsigned s = 0; s < _nStates; ++s ) {
                    _lambda[j]( s, n, 0 ) = ds[s];
                }
            }
        }

        // Converge initial condition entropy variables for One Shot IPM or if truncation order is increased
        if( _settings->GetMaxIterations() == 1 || prevSettings->GetMaxDegree() != _settings->GetMaxDegree() ) {
            prevClosure->SetMaxIterations( 10000 );
            for( unsigned j = 0; j < _nCells; ++j ) {
                // solve on high number of moments if cell is Dirichlet Cell
                if( _mesh->GetBoundaryType( j ) == BoundaryType::DIRICHLET ) {
                    _lambda[j].resize( _nStates, _nMultiElements, _nTotal );
                    _lambda[j].reset( 0, 0, prevSettings->GetNTotal() );
                    _closure->SolveClosureSafe( _lambda[j], u[j], _settings->GetNRefinementLevels() - 1 );
                }
                else
                    prevClosure->SolveClosureSafe( _lambda[j], u[j], prevSettings->GetNRefinementLevels() - 1 );
            }
            prevClosure->SetMaxIterations( maxIterations );
        }
    }

    // for restart with increased number of moments reconstruct solution at finer quad points and compute moments for new truncation order
    if( prevSettings->GetMaxDegree() != _settings->GetMaxDegree() ) {
        MatTens uQFullProc = MatTens( _nCells, Tensor( _nStates, _nMultiElements, _settings->GetNQTotal() ) );
        if( maxIterations == 1 ) _closure->SetMaxIterations( 10000 );    // if one shot IPM is used, make sure that initial duals are converged
        for( unsigned j = 0; j < _nCells; ++j ) {
            // compute solution with low moment order at fine Quadrature nodes
            _closure->U( uQFullProc[j], _closure->EvaluateLambda( _lambda[j] ) );

            // compute higher order moments
            auto uCurrent = uQFullProc[j] * _closure->GetPhiTildeWf();
            u[j].resize( _nStates, _nMultiElements, _nTotal );
            u[j] = uCurrent;    // new Moments of size _nTotal

            // compute lambda with size _nTotal
            if( _lambda[j].columns() != _nTotal ) {
                _lambda[j].resize( _nStates, _nMultiElements, _nTotal );
                _lambda[j].reset( 0, 0, prevSettings->GetNTotal() );
            }
            _closure->SolveClosureSafe( _lambda[j], u[j], _settings->GetNRefinementLevels() - 1 );
        }
        _closure->SetMaxIterations( maxIterations );
    }
    // delete reload closures and settings
    if( prevSettings != _settings ) delete prevSettings;
    if( prevClosure != _closure ) delete prevClosure;
}

MatTens MomentSolver::SetupIC() const {

    MatTens u( _nCells, Tensor( _nStates, _nMultiElements, _nTotal ) );
    Vector xiEta( _settings->GetNDimXi() );
    Matrix uIC( _nStates, _nQTotal );
    Matrix phiTildeWf = _closure->GetPhiTildeWfAtRef( _settings->GetNRefinementLevels() - 1, true );
    std::vector<Vector> IC;
    auto grid   = _closure->GetQuadratureGrid();
    auto xiQuad = grid->GetNodes();
    Vector xiGrid( _nMultiElements, 0.0 );
    double a = -1.0;
    double b = 1.0;
    for( unsigned n = 0; n < _nMultiElements; ++n ) {
        xiGrid[n] = a + n * ( b - a ) / _nMultiElements;
    }
    double dXiGrid = ( b - a ) / _nMultiElements;

    if( _settings->HasICFile() ) {
        IC = _mesh->Import();
    }
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned n = 0; n < _nMultiElements; ++n ) {
            for( unsigned k = 0; k < _nQTotal; ++k ) {
                for( unsigned m = 0; m < _settings->GetNDimXi(); ++m ) {
                    xiEta[m] = xiGrid[n] + 0.5 * ( xiQuad[k][m] + 1.0 ) * dXiGrid;    // TODO: check if not uniform distribution
                }

                if( _settings->HasICFile() ) {
                    column( uIC, k ) = _problem->LoadIC( IC[j], xiEta );
                }
                else {
                    column( uIC, k ) = _problem->IC( _mesh->GetCenterPos( j ), xiEta );
                }
            }

            auto uElement = uIC * phiTildeWf;
            for( unsigned s = 0; s < _nStates; ++s ) {
                for( unsigned i = 0; i < _nTotal; ++i ) {
                    u[j]( s, n, i ) = uElement( s, i );
                }
            }
        }
    }
    return u;
}

void MomentSolver::Export( const MatTens& u, const MatTens& lambda, const VectorU& refinementLevel ) const {
    std::shared_ptr<spdlog::logger> moment_writer = spdlog::get( "moments" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned j = 0; j < _nStates; ++j ) {
                for( unsigned k = 0; k < _nTotal; ++k ) {
                    line << std::setprecision( std::numeric_limits<double>::digits10 ) << u[i]( j, l, k ) << ",";
                }
            }
        }
        moment_writer->info( line.str() );
    }
    moment_writer->flush();
    std::shared_ptr<spdlog::logger> dual_writer = spdlog::get( "duals" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned j = 0; j < _nStates; ++j ) {
                for( unsigned k = 0; k < _nTotal; ++k ) {
                    line << std::setprecision( std::numeric_limits<double>::digits10 ) << lambda[i]( j, l, k ) << ",";
                }
            }
        }
        dual_writer->info( line.str() );
    }
    dual_writer->flush();
    if( _settings->GetProblemType() == P_EULER_2D ) {
        Vector tmp( _nStates, 0.0 );
        std::vector<std::vector<Vector>> uTemperature( _nCells, std::vector( _nMultiElements, Vector( _nTotal, 0.0 ) ) );
        for( unsigned j = 0; j < _nCells; ++j ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                Vector TVec( _settings->GetNQTotalForRef( refinementLevel[j] ) );
                Matrix phiTildeWf = _closure->GetPhiTildeWfAtRef( refinementLevel[j], true );
                for( unsigned k = 0; k < _settings->GetNQTotalForRef( refinementLevel[j] ); ++k ) {
                    _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], l, k, _nTotalForRef[refinementLevel[j]] ) );
                    double rho  = tmp[0];
                    double rhoU = tmp[1];
                    double rhoV = tmp[2];
                    double rhoE = tmp[3];
                    double p    = 0.4 * ( rhoE - 0.5 * ( pow( rhoU, 2 ) + pow( rhoV, 2 ) ) / rho );
                    TVec[k]     = p / ( rho * 287.87 );
                }
                uTemperature[j][l] = phiTildeWf.transpose() * TVec;
            }
        }
        std::shared_ptr<spdlog::logger> T_writer = spdlog::get( "temperature" );
        for( unsigned i = 0; i < _nCells; ++i ) {
            std::stringstream line;
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                for( unsigned k = 0; k < _nTotal; ++k ) {
                    line << std::setprecision( std::numeric_limits<double>::digits10 ) << uTemperature[i][l][k] << ",";
                }
                T_writer->info( line.str() );
            }
        }
        T_writer->flush();
    }
}

void MomentSolver::Export( const MatTens& u, const MatTens& lambda ) const {
    std::shared_ptr<spdlog::logger> moment_writer = spdlog::get( "moments" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned j = 0; j < _nStates; ++j ) {
                for( unsigned k = 0; k < _nTotal; ++k ) {
                    line << std::setprecision( std::numeric_limits<double>::digits10 ) << u[i]( j, l, k ) << ",";
                }
            }
        }
        moment_writer->info( line.str() );
    }
    moment_writer->flush();
    std::shared_ptr<spdlog::logger> dual_writer = spdlog::get( "duals" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned j = 0; j < _nStates; ++j ) {
                for( unsigned k = 0; k < _nTotal; ++k ) {
                    line << std::setprecision( std::numeric_limits<double>::digits10 ) << lambda[i]( j, l, k ) << ",";
                }
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
            // while( line.find( "==================================" ) != std::string::npos ) {
            //    std::getline( file, line );
            //}
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

    return new Settings( inputStream );
}

MatTens MomentSolver::ImportPrevMoments( unsigned nPrevTotal ) const {
    MatTens u( _nCells, Tensor( _nStates, _nMultiElements, nPrevTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_moments" );
    std::string line;
    for( unsigned j = 0; j < _nCells; ++j ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                for( unsigned i = 0; i < nPrevTotal; ++i ) {
                    std::getline( lineStream, cell, ',' );
                    u[j]( s, l, i ) = std::stod( cell );
                }
            }
        }
    }
    file.close();
    return u;
}

MatTens MomentSolver::ImportPrevDuals( unsigned nPrevTotal ) {
    MatTens lambda( _nCells, Tensor( _nStates, _nMultiElements, nPrevTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_duals" );
    std::string line;
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned l = 0; l < _nMultiElements; ++l ) {
                for( unsigned k = 0; k < nPrevTotal; ++k ) {
                    std::getline( lineStream, cell, ',' );
                    lambda[i]( j, l, k ) = std::stod( cell );
                }
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
    if( _settings->GetMeshDimension() == 2 ) {
        for( unsigned j = 0; j < _nCells; ++j ) {
            if( _mesh->GetGrid()[j]->GetCenter()[0] > a[0] && _mesh->GetGrid()[j]->GetCenter()[0] < b[0] &&
                _mesh->GetGrid()[j]->GetCenter()[1] > a[1] && _mesh->GetGrid()[j]->GetCenter()[1] < b[1] ) {

                switch( LNorm ) {
                    case 0:
                        for( unsigned s = 0; s < 2 * _nStates; ++s ) {
                            error[s]   = std::max( std::fabs( ( solution( s, j ) - _referenceSolution[j][s] ) ) * _mesh->GetArea( j ), error[s] );
                            refNorm[s] = std::max( std::fabs( _referenceSolution[j][s] ) * _mesh->GetArea( j ), refNorm[s] );
                        }
                        break;
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
    }
    else {
        for( unsigned j = 0; j < _nCells; ++j ) {
            switch( LNorm ) {
                case 0:
                    for( unsigned s = 0; s < 2 * _nStates; ++s ) {
                        error[s]   = std::max( std::fabs( ( solution( s, j ) - _referenceSolution[j][s] ) ) * _mesh->GetArea( j ), error[s] );
                        refNorm[s] = std::max( std::fabs( _referenceSolution[j][s] ) * _mesh->GetArea( j ), refNorm[s] );
                    }
                    break;
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
        error[s] = std::pow( error[s] / refNorm[s], 1.0 / std::max( double( LNorm ), 1.0 ) );
    }
    return error;
}

void MomentSolver::WriteErrors( const VectorU& refinementLevel, double t ) {
    // define rectangle for error computation
    Vector a( 2 );
    a[0] = -0.05;
    a[1] = -0.5;
    Vector b( 2 );
    b[0]                  = 1.05;
    b[1]                  = 0.5;
    auto l1ErrorMeanLog   = spdlog::get( "l1ErrorMean" );
    auto l2ErrorMeanLog   = spdlog::get( "l2ErrorMean" );
    auto lInfErrorMeanLog = spdlog::get( "lInfErrorMean" );
    auto l1ErrorVarLog    = spdlog::get( "l1ErrorVar" );
    auto l2ErrorVarLog    = spdlog::get( "l2ErrorVar" );
    auto lInfErrorVarLog  = spdlog::get( "lInfErrorVar" );

    Matrix meanAndVar = WriteMeanAndVar( refinementLevel, t, true, false );

    auto l1Error   = this->CalculateError( meanAndVar, 1, a, b );
    auto l2Error   = this->CalculateError( meanAndVar, 2, a, b );
    auto lInfError = this->CalculateError( meanAndVar, 0, a, b );

    std::ostringstream osL1ErrorMean, osL2ErrorMean, osLInfErrorMean, osL1ErrorVar, osL2ErrorVar, osLInfErrorVar;
    for( unsigned i = 0; i < _nStates; ++i ) {
        osL1ErrorMean << std::scientific << l1Error[i] << "\t";
        osL2ErrorMean << std::scientific << l2Error[i] << "\t";
        osLInfErrorMean << std::scientific << lInfError[i] << "\t";
        osL1ErrorVar << std::scientific << l1Error[i + _nStates] << "\t";
        osL2ErrorVar << std::scientific << l2Error[i + _nStates] << "\t";
        osLInfErrorVar << std::scientific << lInfError[i + _nStates] << "\t";
    }

    l1ErrorMeanLog->info( osL1ErrorMean.str() );
    l2ErrorMeanLog->info( osL2ErrorMean.str() );
    lInfErrorMeanLog->info( osLInfErrorMean.str() );
    l1ErrorVarLog->info( osL1ErrorVar.str() );
    l2ErrorVarLog->info( osL2ErrorVar.str() );
    lInfErrorVarLog->info( osLInfErrorVar.str() );
}

void MomentSolver::PerformInitialStep( const VectorU& refinementLevel, MatTens& u ) {
// initial dual solve
#pragma omp parallel for schedule( dynamic, 10 )
    for( unsigned j = 0; j < static_cast<unsigned>( _cellIndexPE.size() ); ++j ) {
        _closure->SolveClosureSafe( _lambda[_cellIndexPE[j]], u[_cellIndexPE[j]], refinementLevel[_cellIndexPE[j]] );
    }

    MatTens uQ = MatTens( _nCells + 1, Tensor( _nStates, _nMultiElements, _settings->GetNqPE() ) );
    // compute solution at quad points
    for( unsigned j = 0; j < _nCells; ++j ) {
        uQ[j] = _closure->U( _closure->EvaluateLambdaOnPE( _lambda[j], refinementLevel[j], refinementLevel[j] ) );
    }

    // compute partial moment vectors on each PE (for inexact dual variables)
    for( unsigned j = 0; j < _nCells; ++j ) {
        u[j] = uQ[j] * _closure->GetPhiTildeWfAtRef( refinementLevel[j] );
    }
}
