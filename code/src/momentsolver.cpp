#include "momentsolver.h"

MomentSolver::MomentSolver( Problem* problem ) : _problem( problem ) {
    _quad = new Legendre( _problem->GetNQuadPoints() );

    _mesh = _problem->GetMesh();

    _x      = _mesh->GetGrid();
    _dx     = _mesh->GetSpacing()[0];    // only equidistant!
    _nCells = _mesh->GetNumCells();
    _a      = _x[0];
    _b      = _x[_x.size() - 1];

    _uL = 12.0;
    _uR = 3.0;

    _nMoments = _problem->GetNMoments();
    _tEnd     = _problem->GetTEnd();
    _nStates  = _problem->GetNStates();

    _closure = Closure::Create( _problem );
    _limiter = Limiter::Create( _closure, _problem );
    _time    = TimeSolver::Create( _problem, _closure );
    _plot    = PlotEngine::Create( _problem );

    _dt = _time->GetTimeStepSize();
}

MomentSolver::~MomentSolver() {
    delete _quad;
    delete _closure;
    delete _limiter;
    delete _time;
    delete _plot;
}

void MomentSolver::Solve() {
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    double t = 0;
    // create solution fields
    std::vector<blaze::DynamicMatrix<double>> uNew( _nCells + 4, blaze::DynamicMatrix<double>( _nStates, _nMoments, 0.0 ) );
    std::vector<blaze::DynamicMatrix<double>> u( _nCells + 4, blaze::DynamicMatrix<double>( _nStates, _nMoments, 0.0 ) );
    _lambda = std::vector<blaze::DynamicMatrix<double>>( _nCells + 4, blaze::DynamicMatrix<double>( _nStates, _nMoments, 0.0 ) );

    u = SetupIC();

    for( int j = 0; j < _nCells + 4; ++j ) {
        _lambda[j] = _closure->SolveClosure( u[j], _lambda[j] );
        u[j]       = CalculateMoments( _lambda[j] );
    }

    // Begin time loop
    while( t < _tEnd ) {
        // Modify moments into realizable direction
        for( int j = 3; j < _nCells + 1; ++j ) {
            u[j] = CalculateMoments( _lambda[j] );
        }

        // Time Update Moments
        _time->Advance(
            std::bind( &MomentSolver::numFlux, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 ),
            uNew,
            u,
            _lambda );

        // Time Update dual variables
        for( int j = 3; j < _nCells + 1; ++j ) {
            _lambda[j] = _closure->SolveClosure( uNew[j], _lambda[j] );
        }
        if( t == 0 )
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
}

blaze::DynamicMatrix<double> MomentSolver::numFlux( const blaze::DynamicMatrix<double>& lambda0,
                                                    const blaze::DynamicMatrix<double>& lambda1,
                                                    const blaze::DynamicMatrix<double>& lambda2,
                                                    const blaze::DynamicMatrix<double>& lambda3 ) {
    blaze::DynamicMatrix<double> g =
        _problem->G( _closure->U( _closure->EvaluateLambda( lambda1 ) ) + 0.5 * _dx * _limiter->Slope( lambda0, lambda1, lambda2 ),
                     _closure->U( _closure->EvaluateLambda( lambda2 ) ) - 0.5 * _dx * _limiter->Slope( lambda1, lambda2, lambda3 ) );

    return 0.5 * g * _closure->GetPhiTildeW();
}

blaze::DynamicMatrix<double> MomentSolver::CalculateMoments( const blaze::DynamicMatrix<double>& lambda ) {
    return 0.5 * _closure->U( _closure->EvaluateLambda( lambda ) ) * _closure->GetPhiTildeW();
}

std::vector<blaze::DynamicMatrix<double>> MomentSolver::SetupIC() {
    std::vector<blaze::DynamicMatrix<double>> out( _nCells + 4, blaze::DynamicMatrix<double>( _nStates, _nMoments, 0.0 ) );
    blaze::DynamicVector<double> xi       = _quad->GetNodes();
    blaze::DynamicVector<double> w        = _quad->GetWeights();
    blaze::DynamicMatrix<double> phiTilde = _closure->GetPhiTilde();
    for( int j = 0; j < _nCells + 4; ++j ) {
        for( int l = 0; l < _nStates; ++l ) {
            for( int i = 0; i < _nMoments; ++i ) {
                for( int k = 0; k < _problem->GetNQuadPoints(); ++k ) {
                    out[j]( l, i ) += 0.5 * w[k] * IC( _x[j], xi[k] ) * phiTilde( k, i );
                }
            }
        }
    }
    return out;
}

double MomentSolver::IC( double x, double xi ) {
    double a     = 0.5;
    double b     = 1.5;
    double sigma = 0.7;
    if( x < a + sigma * xi ) {
        return _uL;
    }
    else if( x < b + sigma * xi ) {
        return _uL + ( _uR - _uL ) * ( a + sigma * xi - x ) / ( a - b );
    }
    else {
        return _uR;
    }
}

void MomentSolver::Plot() {
    int cellIndex, nXi;
    int plotState = 0;
    try {
        auto file = cpptoml::parse_file( _problem->GetInputFile() );
        auto plot = file->get_table( "plot" );

        cellIndex = plot->get_as<int>( "cellIndex" ).value_or( -1 );
        nXi       = plot->get_as<int>( "evalPoints" ).value_or( -1 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _problem->GetInputFile() << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }

    blaze::DynamicVector<double> xi( nXi, 0.0 );
    for( int k = 0; k < nXi; ++k ) {
        xi[k] = -1 + k / ( nXi - 1.0 ) * 2.0;
    }
    blaze::DynamicMatrix<double> resM    = _closure->U( _closure->EvaluateLambda( _lambda[cellIndex - 1], xi ) );
    blaze::DynamicMatrix<double> resQuad = _closure->U( _closure->EvaluateLambda( _lambda[cellIndex - 1] ) );

    // save first state and calculate exact solution
    blaze::DynamicVector<double> res       = blaze::DynamicVector<double>( nXi, 0.0 );
    blaze::DynamicVector<double> resCourse = blaze::DynamicVector<double>( _problem->GetNQuadPoints(), 0.0 );
    blaze::DynamicVector<double> exRes( nXi );
    for( int k = 0; k < nXi; ++k ) {
        res[k]   = resM( plotState, k );
        exRes[k] = _problem->ExactSolution( _tEnd, _x[cellIndex - 1], xi[k] );
    }

    for( int k = 0; k < _problem->GetNQuadPoints(); ++k ) {
        resCourse[k] = resQuad( plotState, k );
    }

    _plot->Plot1D( xi, res, xi, exRes );
    //_plot->Plot1D(_quad->GetNodes(), resCourse, xi, exRes);
}

void MomentSolver::PlotFixedXi() {
    double xi     = 0.0;
    int plotState = 0;
    blaze::DynamicVector<double> xiVec( _nCells + 4, xi );

    int nFine;
    try {
        auto file = cpptoml::parse_file( _problem->GetInputFile() );
        auto plot = file->get_table( "plot" );

        nFine = plot->get_as<int>( "evalPoints" ).value_or( -1 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _problem->GetInputFile() << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }

    blaze::DynamicMatrix<double> resM( _nCells + 4, _nStates );
    blaze::DynamicVector<double> resVec( _nStates, 0.0 );
    for( int j = 0; j < _nCells + 4; ++j ) {
        _closure->U( resVec, _closure->EvaluateLambda( _lambda[j], xiVec, j ) );
        row( resM, j ) = blaze::trans( resVec );
    }

    blaze::DynamicVector<double> res( _nCells + 4 );
    for( int j = 0; j < _nCells + 4; ++j ) {
        res[j] = resM( j, plotState );
    }

    blaze::DynamicVector<double> xFine( nFine, 0.0 ), exRes( nFine, 0.0 );

    for( int j = 0; j < nFine; ++j ) {
        xFine[j] = _a + j * ( _b - _a ) / ( nFine - 1 );
        exRes[j] = _problem->ExactSolution( _tEnd, xFine[j], xi );
    }

    _plot->Plot1D( _x, res, xFine, exRes );
}

void MomentSolver::PlotExpectedValue() {
    const unsigned int nQuadFine            = 200;
    Legendre* quadFine                      = new Legendre( nQuadFine );
    blaze::DynamicVector<double> wFine      = quadFine->GetWeights();
    blaze::DynamicVector<double> xiQuadFine = quadFine->GetNodes();
    blaze::DynamicVector<double> w          = _quad->GetWeights();
    blaze::DynamicVector<double> xiQuad     = _quad->GetNodes();
    unsigned int nFine;

    try {
        auto file = cpptoml::parse_file( _problem->GetInputFile() );
        auto plot = file->get_table( "plot" );

        nFine = plot->get_as<int>( "evalPoints" ).value_or( -1 );
    } catch( const cpptoml::parse_exception& e ) {
        std::cerr << "Failed to parse " << _problem->GetInputFile() << ": " << e.what() << std::endl;
        exit( EXIT_FAILURE );
    }

    blaze::DynamicVector<double> xFine( nFine, 0.0 ), exRes( nFine, 0.0 ), res( _nCells + 4, 0.0 );

    for( unsigned int j = 0; j < nFine; ++j ) {
        xFine[j] = _a + j * ( _b - _a ) / ( nFine - 1 );
        for( unsigned int k = 0; k < nQuadFine; ++k ) {
            exRes[j] += 0.5 * wFine[k] * _problem->ExactSolution( _tEnd, xFine[j], xiQuadFine[k] );
        }
    }

    blaze::DynamicVector<double> resVec( _nStates, 0.0 );
    for( int j = 0; j < _nCells + 4; ++j ) {
        for( int k = 0; k < _problem->GetNQuadPoints(); ++k ) {
            _closure->U( resVec, _closure->EvaluateLambda( _lambda[j], xiQuad, k ) );
            res[j] += 0.5 * w[k] * resVec[0];
        }
    }

    _plot->Plot1D( _x, res, xFine, exRes );
}

void MomentSolver::Print() {
    /*
    int cellIndex = 55;
    int nXi = 20;
    blaze::DynamicVector<double> xi(nXi,0.0);
    for( int k = 0; k<nXi; ++k ){
        xi[k] = -1 + k/(nXi-1.0)*2.0;
    }
    std::cout<<_x[cellIndex]<<std::endl;
    std::cout<<_closure->U(_closure->EvaluateLambda(_lambda[cellIndex-1],xi))<<std::endl;*/
}
