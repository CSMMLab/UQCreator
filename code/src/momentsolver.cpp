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
}

MomentSolver::~MomentSolver() {
    delete _closure;
    delete _time;
}

void MomentSolver::Solve() {
    auto log = spdlog::get( "event" );

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    MatVec uQ = MatVec( _nCells + 1, Matrix( _nStates, _settings->GetNQTotal() ) );

    // create solution fields
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    if( _settings->HasRestartFile() ) {
        this->ImportTime();
        uQ = this->ImportMoments();
    }
    else {
        uQ = SetupIC();
    }
    MatVec uNew = uQ;
    MatVec uOld = uNew;

    Vector ds( _nStates );
    Vector u0( _nStates );

    double residualFull;

    std::vector<unsigned> cellIndexPE = _settings->GetCellIndexPE();
    std::vector<int> PEforCell        = _settings->GetPEforCell();

    auto numFluxPtr = std::bind( &MomentSolver::numFlux,
                                 this,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5 );

    std::cout << uQ[2]( 0, 0 ) << std::endl;

    log->info( "{:10}   {:10}", "t", "residual" );
    // Begin time loop
    double t, dt;
    for( t = _tStart; t < _tEnd; ) {
        double residual = 0;

        // determine time step size
        double dtCurrent;
        dt = 1e10;
        for( unsigned j = 0; j < _nCells; ++j ) {
            dtCurrent = _problem->ComputeDt( uQ[j], _mesh->GetMaxEdge( j ) / _mesh->GetArea( j ) );
            if( dtCurrent < dt ) dt = dtCurrent;
        }
        t += dt;

        _time->Advance( numFluxPtr, uNew, uQ, dt );

        for( unsigned j = 0; j < _nCells; ++j ) {
            residual += std::fabs( uNew[j]( 0, 0 ) - uQ[j]( 0, 0 ) );    // * _mesh->GetArea( j ) / dt;
        }
        log->info( "{:03.8f}   {:01.5e}", t, residual );
        std::cout << "dt = " << dt << std::endl;

        // update solution
        for( unsigned j = 0; j < _nCells; ++j ) {
            uQ[j] = uNew[j];
        }
    }

    // save final moments on u
    uQ = uNew;

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
            for( unsigned i = 0; i < _nStates; ++i ) {
                meanAndVar( i, j ) = uQ[j]( i, 0 );
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
                    meanAndVar( 2 * _nStates + i, j ) += tmp[i];
                }
            }
        }
    }

    if( _settings->GetProblemType() == P_SHALLOWWATER_2D )
        _mesh->ExportShallowWater( meanAndVar );
    else
        _mesh->Export( meanAndVar );

    this->Export( uNew, _lambda );

    // for( unsigned j = 0; j < _nCells; ++j ) std::cout << uQ[j] << std::endl;
}

void MomentSolver::numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    out += _problem->G( u1, u2, nUnit, n );
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

        u[j] = uIC;
    }
    return u;
}

void MomentSolver::Export( const MatVec& u, const MatVec& lambda ) const {
    std::shared_ptr<spdlog::logger> moment_writer = spdlog::get( "moments" );
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::stringstream line;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < _nTotal; ++k ) {
                line << u[i]( j, k ) << ",";
            }
        }
        moment_writer->info( line.str() );
    }
    moment_writer->flush();
}

void MomentSolver::ImportTime() {
    auto file = std::ifstream( _settings->GetRestartFile() );
    std::string line;
    while( std::getline( file, line ) ) {
        if( line.find( "tEnd" ) != std::string::npos ) {
            line.erase( 0, line.find_first_of( '=' ) + 1 );
            line.erase( std::remove_if( line.begin(), line.end(), []( char c ) -> bool { return std::isspace<char>( c, std::locale::classic() ); } ),
                        line.end() );
            _tStart = std::stod( line );
            break;
        }
    }
}

MatVec MomentSolver::ImportMoments() {
    MatVec u( _nCells, Matrix( _nStates, _nTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_moments" );
    std::string line;
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

MatVec MomentSolver::ImportDuals() {
    MatVec lambda( _nCells + 1, Matrix( _nStates, _nTotal ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_duals" );
    std::string line;
    for( unsigned i = 0; i < _nCells + 1; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < _nTotal; ++k ) {
                std::getline( lineStream, cell, ',' );
                lambda[i]( j, k ) = std::stod( cell );
            }
        }
    }
    file.close();
    return lambda;
}
