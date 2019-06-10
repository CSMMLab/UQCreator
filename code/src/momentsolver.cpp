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

MatVec MomentSolver::Solve( Vector xi ) {
    auto log = spdlog::get( "event" );

    MatVec uQ = MatVec( _nCells + 1, Matrix( _nStates, 1 ) );

    // create solution fields
    MatVec u( _nCells, Matrix( _nStates, 1 ) );
    if( _settings->HasRestartFile() ) {
        this->ImportTime();
        uQ = this->ImportMoments();
    }
    else {
        uQ = SetupIC( xi );
    }
    MatVec uNew = uQ;

    auto numFluxPtr = std::bind( &MomentSolver::numFlux,
                                 this,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5 );

    log->info( "{:10}   {:10}", "t", "residual" );
    // Begin time loop
    double t = _tStart;
    double dt;
    double minResidual = _settings->GetMinResidual();
    double residual    = minResidual + 1.0;
    while( t < _tEnd && residual > minResidual ) {
        residual = 0;

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
            residual += std::fabs( uNew[j]( 0, 0 ) - uQ[j]( 0, 0 ) ) * _mesh->GetArea( j );
        }
        log->info( "{:03.8f}   {:01.5e}   {:01.5e}", t, residual, residual / dt );
        // std::cout << "dt = " << dt << std::endl;

        // update solution
        for( unsigned j = 0; j < _nCells; ++j ) {
            uQ[j] = uNew[j];
        }
    }

    // save final moments on uQ
    uQ = uNew;

    return uQ;
}

void MomentSolver::numFlux( Matrix& out, const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    out += _problem->G( u1, u2, nUnit, n );
}

MatVec MomentSolver::SetupIC( const Vector& xi ) {
    MatVec u( _nCells, Matrix( _nStates, 1 ) );
    Matrix uIC( _nStates, _nQTotal );
    std::vector<Vector> IC;
    if( _settings->HasICFile() ) {
        IC = _mesh->Import();
    }
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _settings->HasICFile() ) {
            column( uIC, 0 ) = _problem->LoadIC( IC[j], xi );
        }
        else {
            column( uIC, 0 ) = _problem->IC( _mesh->GetCenterPos( j ), xi );
        }

        u[j] = uIC;
    }
    return u;
}

void MomentSolver::Export( const MatVec& u ) const {
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
    MatVec u( _nCells, Matrix( _nStates, 1 ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_moments" );
    std::string line;
    for( unsigned i = 0; i < _nCells; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < 1; ++k ) {
                std::getline( lineStream, cell, ',' );
                u[i]( j, 0 ) = std::stod( cell );
            }
        }
    }
    file.close();
    return u;
}

MatVec MomentSolver::ImportDuals() {
    MatVec lambda( _nCells + 1, Matrix( _nStates, 1 ) );
    auto file = std::ifstream( _settings->GetRestartFile() + "_duals" );
    std::string line;
    for( unsigned i = 0; i < _nCells + 1; ++i ) {
        std::getline( file, line );
        std::stringstream lineStream( line );
        std::string cell;
        for( unsigned j = 0; j < _nStates; ++j ) {
            for( unsigned k = 0; k < 1; ++k ) {
                std::getline( lineStream, cell, ',' );
                lambda[i]( j, 0 ) = std::stod( cell );
            }
        }
    }
    file.close();
    return lambda;
}
