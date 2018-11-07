#include "momentsolver.h"

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
}

MomentSolver::~MomentSolver() {
    delete _quad;
    delete _closure;
    delete _time;
}

void MomentSolver::Solve() {
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    // create solution fields
    MatVec u    = SetupIC();
    MatVec uNew = MatVec( _nCells, Matrix( _nStates, _nMoments ) );
    MatVec uQ   = MatVec( _nCells + 1, Matrix( _nStates, _nQuadPoints ) );
    _lambda     = _problem->InitLambda( u );

    auto numFluxPtr = std::bind( &MomentSolver::numFlux,
                                 this,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::placeholders::_4,
                                 std::placeholders::_5 );

#pragma omp parallel for
    for( unsigned j = 0; j < _nCells; ++j ) {
        _closure->SolveClosure( _lambda[j], u[j] );
    }

    std::cout << std::fixed << std::setprecision( 8 ) << std::setw( 10 ) << "t" << std::setw( 15 ) << "residual" << std::endl;
    // Begin time loop
    for( double t = 0.0; t < _tEnd; t += _dt ) {
        double residual = 0;
        if( t == 0.0 ) {
#pragma omp parallel for
            for( unsigned j = 0; j < _nCells; ++j ) {
                uQ[j] = _closure->U( _closure->EvaluateLambda( _lambda[j] ) );
                u[j]  = 0.5 * uQ[j] * _closure->GetPhiTildeW();
            }
        }
        else {
#pragma omp parallel for reduction( + : residual )
            for( unsigned j = 0; j < _nCells; ++j ) {
                _closure->SolveClosure( _lambda[j], uNew[j] );
                residual += std::fabs( uNew[j]( 0, 0 ) - u[j]( 0, 0 ) ) * _mesh->GetArea( j ) / _dt;

                uQ[j] = _closure->U( _closure->EvaluateLambda( _lambda[j] ) );
                u[j]  = 0.5 * uQ[j] * _closure->GetPhiTildeW();
            }
            std::cout << std::fixed << std::setprecision( 8 ) << std::setw( 10 ) << t << std::setw( 15 ) << residual << std::endl;
        }

        _time->Advance( numFluxPtr, uNew, u, uQ );

        // exit( EXIT_FAILURE );
    }

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
                meanAndVar( i, j ) += 0.5 * w[k] * tmp[i];
            }
        }

        // var
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            _closure->U( tmp, _closure->EvaluateLambda( _lambda[j], k ) );
            for( unsigned i = 0; i < _nStates; ++i ) {
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

MatVec MomentSolver::SetupIC() {
    MatVec u( _nCells, Matrix( _nStates, _nMoments, 0.0 ) );
    Vector xi = _quad->GetNodes();
    Matrix uIC( _nStates, _nQuadPoints, 0.0 );
    Matrix phiTildeW = _closure->GetPhiTildeW();
    for( unsigned j = 0; j < _nCells; ++j ) {
        for( unsigned k = 0; k < _nQuadPoints; ++k ) {
            column( uIC, k ) = IC( _mesh->GetCenterPos( j ), xi[k] );
        }
        u[j] = 0.5 * uIC * phiTildeW;
    }
    return u;
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
        double sigma = 1.25;
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
