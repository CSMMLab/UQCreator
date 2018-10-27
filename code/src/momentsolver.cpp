#include "momentsolver.h"

MomentSolver::MomentSolver( Settings* settings, Mesh* mesh, Problem* problem ) : _settings( settings ), _mesh( mesh ), _problem( problem ) {
    _nCells      = _settings->GetNumCells();
    _nMoments    = _settings->GetNMoments();
    _tEnd        = _settings->GetTEnd();
    _nStates     = _settings->GetNStates();
    _nQuadPoints = _settings->GetNQuadPoints();

    _quad    = new Legendre( _nQuadPoints );
    _closure = Closure::Create( _settings );
    _limiter = Limiter::Create( _settings, _mesh, _closure );
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
    std::vector<Matrix> uQ = std::vector<Matrix>( _nCells + 1, Matrix( _nStates, _nQuadPoints, 0.0 ) );

    for( unsigned j = 0; j < _nCells; ++j ) {
        if( _settings->GetProblemType() == ProblemType::P_EULER_1D || _settings->GetProblemType() == ProblemType::P_EULER_2D ) {

            double gamma       = _settings->GetGamma();
            double rho         = u[j]( 0, 0 );
            double rhoU        = u[j]( 1, 0 );
            double rhoV        = u[j]( 2, 0 );
            double rhoU2       = rhoU * rhoU / rho;
            double rhoV2       = rhoV * rhoV / rho;
            double rhoE        = u[j]( 3, 0 );
            _lambda[j]( 0, 0 ) = ( rhoU2 + rhoV2 + gamma * ( 2 * rho * rhoE - rhoU2 - rhoV2 ) ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) -
                                 std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );
            _lambda[j]( 1, 0 ) = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
            _lambda[j]( 2, 0 ) = -( ( 2 * rho * rhoV ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );

            _lambda[j]( _settings->GetNStates() - 1, 0 ) = -( rho / ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );
        }
        _lambda[j] = _closure->SolveClosure( u[j], _lambda[j] );
        u[j]       = CalculateMoments( _lambda[j] );    // kann raus!
    }
    std::cout << "Dual variables IC computed." << std::endl;

    // Begin time loop
    unsigned nSteps = 0;
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

        double residual = 0;
        for( unsigned j = 0; j < _nCells; ++j ) {
            residual += std::fabs( uNew[j]( 0, 0 ) - u[j]( 0, 0 ) );
        }
        std::cout << " -> E[rho] residual is " << residual << std::endl;

        // Time Update dual variables
        //#pragma omp parallel for
        for( unsigned j = 0; j < _nCells; ++j ) {
            _lambda[j] = _closure->SolveClosure( uNew[j], _lambda[j] );
        }

        Plot( t, nSteps );

        if( std::fabs( t ) <= std::numeric_limits<double>::epsilon() )
            std::cout << std::fixed << std::setprecision( 8 ) << "t = " << t << std::flush;
        else
            std::cout << std::fixed << std::setprecision( 8 ) << "\r"
                      << "t = " << t << std::flush;
        t += _dt;
        nSteps++;
    }

    double maxRho = 0.0;
    double maxMx  = 0.0;
    double maxMy  = 0.0;
    double maxE   = 0.0;
    for( unsigned j = 0; j < _nCells; ++j ) {
        if( uNew[j]( 0, 0 ) > maxRho ) {
            maxRho = uNew[j]( 0, 0 );
        }
        if( uNew[j]( 0, 0 ) > maxMx ) {
            maxMx = uNew[j]( 1, 0 );
        }
        if( uNew[j]( 0, 0 ) > maxMy ) {
            maxMy = uNew[j]( 2, 0 );
        }
        if( uNew[j]( 0, 0 ) > maxE ) {
            maxE = uNew[j]( 3, 0 );
        }
    }
    std::cout << "Max Rho = " << maxRho << ", Max Mx = " << maxMx << ", Max My = " << maxMy << ", Max maxE " << maxE << std::endl;

    _tEnd = t;

    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    std::cout << "\nFinished!\nRuntime: " << std::setprecision( 3 )
              << std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 << "s" << std::endl;

    Matrix meanAndVar( 2 * _nStates, _mesh->GetNumCells(), 0.0 );
    Vector tmp( _nStates, 0.0 );
    auto xiQuad = _quad->GetNodes();
    Vector w    = _quad->GetWeights();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                // mean
                // std::cout << "tmp = " << tmp[i] << std::endl;
                // std::cout << "w = " << w[k] << std::endl;
                meanAndVar( i, j ) += 0.5 * w[k] * tmp[i];
            }
        }

        // var
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
                // std::cout << "tmp = " << tmp[i] << std::endl;
                // std::cout << "mean = " << meanAndVar( i, j ) << std::endl;
                meanAndVar( i + _nStates, j ) += 0.5 * w[k] * pow( tmp[i] - meanAndVar( i, j ), 2 );
            }
        }
    }

    _mesh->Export( meanAndVar );
}

Matrix MomentSolver::numFlux( const Matrix& u1, const Matrix& u2, const Vector& nUnit, const Vector& n ) {
    Matrix g = _problem->G( u1, u2, nUnit, n );
    return 0.5 * g * _closure->GetPhiTildeW();
}

Matrix MomentSolver::CalculateMoments( const Matrix& lambda ) {
    return 0.5 * _closure->U( _closure->EvaluateLambda( lambda ) ) * _closure->GetPhiTildeW();
}

std::vector<Matrix> MomentSolver::SetupIC() {
    std::vector<Matrix> out( _nCells + 1, Matrix( _nStates, _nMoments, 0.0 ) );
    Vector xi = _quad->GetNodes();
    Matrix uIC( _nStates, _nQuadPoints, 0.0 );
    Matrix phiTildeW = _closure->GetPhiTildeW();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            column( uIC, k ) = IC( _mesh->GetCenterPos( j ), xi[k] );
        }
        out[j] = 0.5 * uIC * phiTildeW;
    }
    return out;
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
        double sigma = 0.5;
        double gamma = 1.4;
        double R     = 287.87;
        double T     = 273.15;
        double p     = 101325.0;
        double Ma    = 0.8;
        double a     = sqrt( gamma * R * T );
        double pi    = 3.14159265359;

        double uMax  = Ma * a;
        double angle = ( 1.25 + sigma * xi ) * ( 2.0 * pi ) / 360.0;
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
