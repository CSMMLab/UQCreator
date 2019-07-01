#include "closure.h"
#include "boundedbarrier.h"
#include "eulerclosure.h"
#include "eulerclosure2d.h"
#include "l2filter.h"
#include "lassofilter.h"
#include "logsin.h"
#include "mathtools.h"
#include "shallowwaterclosure.h"
#include "shallowwaterclosure2d.h"
#include "stochasticgalerkin.h"
#include "tensorizedquadrature.h"
#include "uniformsparsegrid.h"

Closure::Closure( Settings* settings )
    : _settings( settings ), _nMoments( _settings->GetNMoments() ), _nQuadPoints( _settings->GetNQuadPoints() ), _nStates( _settings->GetNStates() ),
      _maxIterations( _settings->GetMaxIterations() ) {
    _log = spdlog::get( "event" );
    // initialize classes: Basis Functions and Quadrature rules are defined for Legendre and Hermite
    _basis.resize( 2 );
    _quad.resize( 2 );
    _basis[0] = Polynomial::Create( _settings, _nMoments, DistributionType::D_LEGENDRE );
    _basis[1] = Polynomial::Create( _settings, _nMoments, DistributionType::D_HERMITE );
    _quad[0]  = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_LEGENDRE );
    _quad[1]  = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_HERMITE );

    // get number of uncertain dimensions
    _numDimXi = _settings->GetNDimXi();

    // get index vector for basis function calculation
    std::vector<std::vector<unsigned>> indices = _settings->GetPolyIndices();
    _nTotal                                    = _settings->GetNTotal();

    // define quadrature
    _wGrid.resize( _settings->GetNRefinementLevels() );
    _nQTotalForRef.resize( _settings->GetNRefinementLevels() );
    _nTotalForRef      = _settings->GetNTotalRefinementLevel();
    VectorU quadLevel  = _settings->GetQuadLevel();
    unsigned oldQLevel = 0;
    _quadGrid          = nullptr;
    for( unsigned rlevel = 0; rlevel < _settings->GetNRefinementLevels(); ++rlevel ) {
        if( oldQLevel != quadLevel[rlevel] ) {    // new quadrature weights needed
            if( _quadGrid ) delete _quadGrid;
            _quadGrid = QuadratureGrid::Create( _settings, quadLevel[rlevel] );
        }
        _wGrid[rlevel]         = _quadGrid->GetWeights();
        _nQTotalForRef[rlevel] = _wGrid[rlevel].size();
        oldQLevel              = quadLevel[rlevel];
    }
    _xiGrid = _quadGrid->GetNodes();
    // give number of quad points per level to settings. Quadrature variables are computed here.
    _settings->SetNQTotalForRef( _nQTotalForRef );

    // set total number of quadrature points
    _nQTotal = _quadGrid->GetNodeCount();

    // compute basis functions evaluated at the quadrature points
    _phiTilde    = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeF   = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeWf  = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeVec = std::vector<Vector>( _nQTotal, Vector( _nTotal, 0.0 ) );

    unsigned n = 0;
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned l = 0; l < _numDimXi; ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                _phiTilde( k, i ) *=
                    _basis[n]->Evaluate( indices[i][l], _xiGrid[k][l] ) / _basis[n]->L2Norm( indices[i][l] );    // sqrt( 2.0 * i + 1.0 );
                _phiTildeF( k, i ) *=
                    _basis[n]->Evaluate( indices[i][l], _xiGrid[k][l] ) / _basis[n]->L2Norm( indices[i][l] ) * _basis[n]->fXi( _xiGrid[k][l] );
            }
            _phiTildeWf( k, i ) = _phiTildeF( k, i ) * _wGrid[_settings->GetNRefinementLevels() - 1][k];
            _phiTildeVec[k][i]  = _phiTilde( k, i );
        }
    }

    _phiTildeTrans      = trans( _phiTilde );
    auto phiTildeFTrans = trans( _phiTildeF );
    // calculate partial matrix for Hessian calculation
    _hPartial = MatVec( _nQTotal, Matrix( _nTotal, _nTotal ) );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        _hPartial[k] = outer( column( _phiTildeTrans, k ), column( phiTildeFTrans, k ) );    // TODO
    }
    /*
        // Test if polynomials are orthonormal
        Matrix testQuad( _nTotal, _nTotal, 0.0 );
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned j = 0; j < _nTotal; ++j ) {
                for( unsigned k = 0; k < _nQTotal; ++k ) {
                    testQuad( i, j ) += _phiTilde( k, i ) * _phiTildeWf( k, j );
                }
            }
        }
        std::cout << "test I " << testQuad << std::endl;*/
}

Closure::~Closure() {
    for( unsigned l = 0; l < _basis.size(); ++l ) {
        delete _basis[l];
        delete _quad[l];
    }
    delete _quadGrid;
}

Closure* Closure::Create( Settings* settings ) {
    auto log         = spdlog::get( "event" );
    auto closureType = settings->GetClosureType();
    if( closureType == ClosureType::C_BOUNDEDBARRIER ) {
        return new BoundedBarrier( settings );
    }
    else if( closureType == ClosureType::C_LOGSIN ) {
        return new LogSin( settings );
    }
    else if( closureType == ClosureType::C_STOCHASTICGALERKIN ) {
        return new StochasticGalerkin( settings );
    }
    else if( closureType == ClosureType::C_EULER_1D ) {
        return new EulerClosure( settings );
    }
    else if( closureType == ClosureType::C_EULER_2D ) {
        return new EulerClosure2D( settings );
    }
    else if( closureType == ClosureType::C_SHALLOWWATER_1D ) {
        return new ShallowWaterClosure( settings );
    }
    else if( closureType == ClosureType::C_SHALLOWWATER_2D ) {
        return new ShallowWaterClosure2D( settings );
    }
    else if( closureType == ClosureType::C_L2FILTER ) {
        return new L2Filter( settings );
    }
    else if( closureType == ClosureType::C_LASSOFILTER ) {
        return new LassoFilter( settings );
    }
    else {
        log->error( "[closure]: Invalid closure type" );
        exit( EXIT_FAILURE );
    }
}

void Closure::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    int maxRefinements = 1000;
    unsigned nTotal    = _nTotalForRef[refLevel];

    Vector g( _nStates * nTotal );

    // check if initial guess is good enough
    Gradient( g, lambda, u, refLevel );
    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
        return;
    }
    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda, refLevel );
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
        return;
    }
    Matrix lambdaNew( _nStates, nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew, nTotal );
    Gradient( dlambdaNew, lambdaNew, u, refLevel );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, u, refLevel );
            dlambda = -g;
            Hessian( H, lambda, refLevel );
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, refLevel );
        }
        int refinementCounter = 0;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, refLevel );
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++refinementCounter > maxRefinements ) {
                _log->error( "[closure] Newton needed too many refinement steps!" );
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

double Closure::CalcNorm( Vector& test, unsigned nTotal ) const {
    double out = 0.0;
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out += pow( test[l * nTotal + i], 2 );
        }
    }
    return sqrt( out );
}

Vector Closure::EvaluateLambda( const Matrix& lambda, unsigned k, unsigned nTotal ) {
    Vector out( _nStates, 0.0 );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out[s] += lambda( s, i ) * _phiTildeVec[k][i];
        }
    }
    return out;
}

Matrix Closure::EvaluateLambda( const Matrix& lambda ) const { return lambda * _phiTildeTrans; }

Matrix Closure::EvaluateLambdaOnPE( const Matrix& lambda, unsigned levelOld, unsigned levelNew ) const {
    Matrix out( _settings->GetNStates(), _settings->GetNqPEAtRef( levelNew ), 0.0 );
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( levelNew );
    unsigned nTotal              = _nTotalForRef[levelOld];

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned k = 0; k < qIndex.size(); ++k ) {
            for( unsigned i = 0; i < nTotal; ++i ) {
                out( s, k ) += lambda( s, i ) * _phiTildeTrans( i, qIndex[k] );
            }
        }
    }
    return out;
}

void Closure::EvaluateLambda( Matrix& out, const Matrix& lambda ) const { out = lambda * _phiTildeTrans; }

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    // std::cout << "Start Hessian" << std::endl;
    Vector uKinetic( _nStates, 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    g.reset();

    // std::cout << "Lambda = " << lambda * _phiTildeVec[nQTotal - 1] << std::endl;

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        // std::cout << "uKinetic = " << uKinetic << std::endl;
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeF( k, i ) * _wGrid[refLevel][k];
            }
        }
    }
    // std::cout << "uKinetic = " << uKinetic << std::endl;
    // std::cout << "g int = " << g << std::endl;

    SubstractVectorMatrixOnVector( g, u, _nTotalForRef[refLevel] );
    // std::cout << "End Gradient" << std::endl;
}

void Closure::Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel ) {
    // std::cout << "Start Hessian" << std::endl;
    // std::cout << "refLevel " << refLevel << ", nQ = " << _nQTotalForRef[refLevel] << ", h size " << _hPartial.size()
    //          << ", nTotal = " << _nTotalForRef[refLevel] << ", size H " << H.columns() << " " << H.rows() << std::endl;
    H.reset();
    Matrix dUdLambda( _nStates, _nStates );    // TODO: preallocate Matrix for Hessian computation -> problems omp
    unsigned nTotal = _nTotalForRef[refLevel];

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dUdLambda, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned m = 0; m < _nStates; ++m ) {
                for( unsigned j = 0; j < nTotal; ++j ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        H( m * nTotal + j, l * nTotal + i ) += _hPartial[k]( j, i ) * _wGrid[refLevel][k] * dUdLambda( l, m );
                    }
                }
            }
        }
    }
    // std::cout << "End Hessian" << std::endl;
}

void Closure::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            y( l, j ) = A( l, j ) + b[l * nTotal + j];
        }
    }
}

void Closure::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            b[l * nTotal + j] = b[l * nTotal + j] - A( l, j );
        }
    }
}

void Closure::DS( Vector& ds, const Vector& u ) const {}

std::vector<Polynomial*> Closure::GetBasis() { return _basis; }

std::vector<Polynomial*> Closure::GetQuadrature() { return _quad; }

void Closure::SetAlpha( double alpha ) { _alpha = alpha; }

void Closure::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    int maxRefinements = 1000;
    unsigned nTotal    = _nTotalForRef[refLevel];

    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector g( _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );

    Vector dlambda = -g;
    Matrix lambdaNew( _nStates, nTotal );

    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        Gradient( g, lambda, u, refLevel );
        dlambda = -g;
        Hessian( H, lambda, refLevel );
        posv( H, g );
        AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
        Gradient( dlambdaNew, lambdaNew, u, refLevel );
        int refinementCounter = 0;
        // std::cout << "Res " << CalcNorm( dlambdaNew, nTotal ) << std::endl;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, refLevel );
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++refinementCounter > maxRefinements ) {
                _log->error( "[closure] Newton needed too many refinement steps!" );
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

void Closure::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }

unsigned Closure::GetMaxIterations() const { return _maxIterations; }

QuadratureGrid* Closure::GetQuadratureGrid() { return _quadGrid; }

Matrix Closure::GetPhiTildeWfAtRef( unsigned level ) const {
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( level );
    Matrix phiTildeWfTrans( unsigned( qIndex.size() ), _nTotalForRef[level], false );
    for( unsigned k = 0; k < qIndex.size(); ++k ) {
        for( unsigned i = 0; i < _nTotalForRef[level]; ++i ) {
            phiTildeWfTrans( k, i ) = _phiTildeF( qIndex[k], i ) * _wGrid[level][qIndex[k]];
        }
    }
    return phiTildeWfTrans;
}

Matrix Closure::GetPhiTildeWfAtRef( unsigned level, bool full ) const {
    unsigned kStart = 0;
    unsigned kEnd   = _nQTotalForRef[level] - 1;
    Matrix phiTildeWfTrans( _nQTotalForRef[level], _nTotalForRef[level], false );
    for( unsigned k = kStart; k <= kEnd; ++k ) {
        for( unsigned i = 0; i < _nTotalForRef[level]; ++i ) {
            phiTildeWfTrans( k - kStart, i ) = _phiTildeF( k, i ) * _wGrid[level][k];
        }
    }
    return phiTildeWfTrans;
}
