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

    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();

    // setup map from k (0,...,_nQTotal-1) to individual indices
    std::vector<std::vector<unsigned>> indices;

    // setup map from i (0,...,nTotal-1) to individual indices
    std::vector<unsigned> indexTest;
    indexTest.resize( _numDimXi );
    unsigned totalDegree;    // compute total degree of basis function i
    for( unsigned i = 0; i < std::pow( _settings->GetNMoments(), _numDimXi ); ++i ) {
        totalDegree = 0;
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indexTest[l] = unsigned( ( i - i % unsigned( std::pow( _nMoments, l ) ) ) / unsigned( std::pow( _nMoments, l ) ) ) % _nMoments;
            totalDegree += indexTest[l];
        }
        // if total degree is sufficiently small or max degree is used, indices are stored
        if( totalDegree < _nMoments || _settings->UsesMaxDegree() ) indices.push_back( indexTest );
    }
    _nTotal = unsigned( indices.size() );

    // define quadrature
    _quadGrid = QuadratureGrid::Create( _settings );
    // TensorizedQuadrature* quadGrid = new TensorizedQuadrature( _settings );
    std::cout << "quadGrid done" << std::endl;
    _xiGrid    = _quadGrid->GetNodes();
    auto wGrid = _quadGrid->GetWeights();

    // set total number of quadrature points
    _nQTotal = _quadGrid->GetNodeCount();

    // compute basis functions evaluated at the quadrature points
    _phiTilde    = Matrix( _nQTotal, _nTotal, 1.0 );
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
                _phiTildeWf( k, i ) *=
                    _basis[n]->Evaluate( indices[i][l], _xiGrid[k][l] ) / _basis[n]->L2Norm( indices[i][l] ) * _basis[n]->fXi( _xiGrid[k][l] );
            }

            _phiTildeWf( k, i ) *= wGrid[k];
            _phiTildeVec[k][i] = _phiTilde( k, i );
        }
    }

    _phiTildeTrans       = trans( _phiTilde );
    auto phiTildeWfTrans = trans( _phiTildeWf );
    // calculate partial matrix for Hessian calculation
    _hPartial = MatVec( _nQTotal, Matrix( _nTotal, _nTotal ) );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        _hPartial[k] = outer( column( _phiTildeTrans, k ), column( phiTildeWfTrans, k ) );    // TODO
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

void Closure::SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal ) {
    int maxRefinements = 1000;

    Vector g( _nStates * nTotal );

    // check if initial guess is good enough
    Gradient( g, lambda, u, nTotal, nQTotal );
    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
        return;
    }
    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );
    // std::cout << "before first Hessian inversion..." << std::endl;
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda, nTotal, nQTotal );
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
        return;
    }
    Matrix lambdaNew( _nStates, nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew, nTotal );
    Gradient( dlambdaNew, lambdaNew, u, nTotal, nQTotal );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, u, nTotal, nQTotal );
            dlambda = -g;
            Hessian( H, lambda, nTotal, nQTotal );
            // std::cout << H << std::endl;
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, nTotal, nQTotal );
        }
        int refinementCounter = 0;
        // std::cout << "Res is " << CalcNorm( dlambda, nTotal ) << std::endl;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, nTotal, nQTotal );
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

Matrix Closure::EvaluateLambdaOnPE( const Matrix& lambda, unsigned nTotal ) const {
    Matrix out( _settings->GetNStates(), _settings->GetNqPE(), 0.0 );
    unsigned kStart = _settings->GetKStart();
    unsigned kEnd   = _settings->GetKEnd();

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned k = kStart; k <= kEnd; ++k ) {
            for( unsigned i = 0; i < nTotal; ++i ) {
                out( s, k - kStart ) += lambda( s, i ) * _phiTildeTrans( i, k );
            }
        }
    }
    return out;
}

void Closure::EvaluateLambda( Matrix& out, const Matrix& lambda ) const { out = lambda * _phiTildeTrans; }

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal ) {
    Vector uKinetic( _nStates, 0.0 );
    g.reset();

    // std::cout << "Lambda = " << lambda * _phiTildeVec[nQTotal - 1] << std::endl;

    for( unsigned k = 0; k < nQTotal; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        // std::cout << "uKinetic = " << uKinetic << std::endl;
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeWf( k, i );
            }
        }
    }
    // std::cout << "uKinetic = " << uKinetic << std::endl;
    // std::cout << "g int = " << g << std::endl;

    SubstractVectorMatrixOnVector( g, u, nTotal );
}

void Closure::Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal ) {
    H.reset();
    Matrix dUdLambda( _nStates, _nStates );    // TODO: preallocate Matrix for Hessian computation -> problems omp

    for( unsigned k = 0; k < nQTotal; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dUdLambda, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned m = 0; m < _nStates; ++m ) {
                for( unsigned j = 0; j < nTotal; ++j ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        H( m * nTotal + j, l * nTotal + i ) += _hPartial[k]( j, i ) * dUdLambda( l, m );
                    }
                }
            }
        }
    }
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

void Closure::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal ) {
    int maxRefinements = 1000;

    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector g( _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );

    Vector dlambda = -g;
    Matrix lambdaNew( _nStates, nTotal );

    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        Gradient( g, lambda, u, nTotal, nQTotal );
        dlambda = -g;
        Hessian( H, lambda, nTotal, nQTotal );
        posv( H, g );
        AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
        Gradient( dlambdaNew, lambdaNew, u, nTotal, nQTotal );
        int refinementCounter = 0;
        // std::cout << "Res " << CalcNorm( dlambdaNew, nTotal ) << std::endl;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, nTotal, nQTotal );
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
