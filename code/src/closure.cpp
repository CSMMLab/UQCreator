#include "closure.h"
#include "boundedbarrier.h"
#include "eulerclosure.h"
#include "eulerclosure2d.h"
#include "l2filter.h"
#include "lassofilter.h"
#include "mathtools.h"
#include "shallowwaterclosure.h"
#include "shallowwaterclosure2d.h"
#include "stochasticgalerkin.h"

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
    _nQTotal  = unsigned( std::pow( _settings->GetNQuadPoints(), _numDimXi ) );

    // setup map from k (0,...,_nQTotal-1) to individual indices
    std::vector<std::vector<unsigned>> indices;
    std::vector<std::vector<unsigned>> indicesQ;
    indicesQ.resize( _nQTotal );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        indicesQ[k].resize( _numDimXi );
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indicesQ[k][l] = unsigned( ( k - k % unsigned( std::pow( _nQuadPoints, l ) ) ) / unsigned( std::pow( _nQuadPoints, l ) ) ) % _nQuadPoints;
        }
    }

    // setup map from i (0,...,_nTotal-1) to individual indices
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

    // store quad points and weights for individual distributions
    std::vector<Vector> xi;
    xi.resize( _quad.size() );
    std::vector<Vector> w;
    w.resize( _quad.size() );
    for( unsigned l = 0; l < _quad.size(); ++l ) {
        xi[l] = _quad[l]->GetNodes();
        w[l]  = _quad[l]->GetWeights();
    }

    // compute basis functions evaluated at the quadrature points
    _phiTilde    = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeWf  = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeVec = std::vector<Vector>( _nQTotal, Vector( _nTotal, 0.0 ) );

    unsigned n;
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned l = 0; l < _numDimXi; ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                _phiTilde( k, i ) *=
                    _basis[n]->Evaluate( indices[i][l], xi[n][indicesQ[k][l]] ) / _basis[n]->L2Norm( indices[i][l] );    // sqrt( 2.0 * i + 1.0 );
                _phiTildeWf( k, i ) *= _basis[n]->Evaluate( indices[i][l], xi[n][indicesQ[k][l]] ) / _basis[n]->L2Norm( indices[i][l] ) *
                                       w[n][indicesQ[k][l]] * _basis[n]->fXi( xi[n][indicesQ[k][l]] );
            }
            // multiplied by pdf
            _phiTildeVec[k][i] = _phiTilde( k, i );    // sqrt( 2.0 * i + 1.0 );
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
        std::cout << testQuad << std::endl;*/
}

Closure::~Closure() {
    for( unsigned l; l < _basis.size(); ++l ) {
        delete _basis[l];
        delete _quad[l];
    }
}

Closure* Closure::Create( Settings* settings ) {
    auto log         = spdlog::get( "event" );
    auto closureType = settings->GetClosureType();
    if( closureType == ClosureType::C_BOUNDEDBARRIER ) {
        return new BoundedBarrier( settings );
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

void Closure::SolveClosure( Matrix& lambda, const Matrix& u ) {
    int maxRefinements = 1000;

    Vector g( _nStates * _nTotal );

    // check if initial guess is good enough
    Gradient( g, lambda, u );
    if( CalcNorm( g ) < _settings->GetEpsilon() ) {
        return;
    }
    Matrix H( _nStates * _nTotal, _nStates * _nTotal );
    Vector dlambdaNew( _nStates * _nTotal );
    // std::cout << "before first Hessian inversion..." << std::endl;
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda );
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda );
        return;
    }
    Matrix lambdaNew( _nStates, _nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew );
    Gradient( dlambdaNew, lambdaNew, u );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, u );
            dlambda = -g;
            Hessian( H, lambda );
            // std::cout << H << std::endl;
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew );
            Gradient( dlambdaNew, lambdaNew, u );
        }
        int refinementCounter = 0;
        // std::cout << "Res is " << CalcNorm( dlambda ) << std::endl;
        while( CalcNorm( dlambda ) < CalcNorm( dlambdaNew ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew );
            Gradient( dlambdaNew, lambdaNew, u );
            if( CalcNorm( dlambdaNew ) < _settings->GetEpsilon() ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++refinementCounter > maxRefinements ) {
                _log->error( "[closure] Newton needed too many refinement steps!" );
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( CalcNorm( dlambdaNew ) < _settings->GetEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

double Closure::CalcNorm( Vector& test ) {
    double out = 0.0;
    double tmp;
    for( unsigned l = 0; l < _nStates; ++l ) {
        tmp = 0.0;
        for( unsigned i = 0; i < _nTotal; ++i ) {
            tmp += pow( test[l * _nTotal + i], 2 );
        }
        out += sqrt( tmp );
    }
    return out;
}

Vector Closure::EvaluateLambda( const Matrix& lambda, unsigned k ) { return lambda * _phiTildeVec[k]; }

Matrix Closure::EvaluateLambda( const Matrix& lambda ) const { return lambda * _phiTildeTrans; }

Matrix Closure::EvaluateLambdaOnPE( const Matrix& lambda ) const {
    Matrix out( _settings->GetNStates(), _settings->GetNqPE(), 0.0 );
    unsigned kStart = _settings->GetKStart();
    unsigned kEnd   = _settings->GetKEnd();

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned k = kStart; k <= kEnd; ++k ) {
            for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
                out( s, k - kStart ) += lambda( s, i ) * _phiTildeTrans( i, k );
            }
        }
    }
    return out;
}

void Closure::EvaluateLambda( Matrix& out, const Matrix& lambda ) const { out = lambda * _phiTildeTrans; }

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u ) {
    Vector uKinetic( _nStates, 0.0 );
    g.reset();

    // std::cout << "Lambda = " << lambda * _phiTildeVec[_nQTotal - 1] << std::endl;

    for( unsigned k = 0; k < _nQTotal; ++k ) {
        U( uKinetic, lambda * _phiTildeVec[k] );
        // std::cout << "uKinetic = " << uKinetic << std::endl;
        for( unsigned j = 0; j < _nTotal; ++j ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * _nTotal + j] += uKinetic[l] * _phiTildeWf( k, j );
            }
        }
    }
    // std::cout << "uKinetic = " << uKinetic << std::endl;
    // std::cout << "g int = " << g << std::endl;

    SubstractVectorMatrixOnVector( g, u );
}

void Closure::Hessian( Matrix& H, const Matrix& lambda ) {
    H.reset();
    Matrix dUdLambda( _nStates, _nStates );    // TODO: preallocate Matrix for Hessian computation -> problems omp

    for( unsigned k = 0; k < _nQTotal; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dUdLambda, lambda * _phiTildeVec[k] );
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned m = 0; m < _nStates; ++m ) {
                for( unsigned j = 0; j < _nTotal; ++j ) {
                    for( unsigned i = 0; i < _nTotal; ++i ) {
                        H( m * _nTotal + j, l * _nTotal + i ) += _hPartial[k]( j, i ) * dUdLambda( l, m );
                    }
                }
            }
        }
    }
}

void Closure::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < _nTotal; ++j ) {
            y( l, j ) = A( l, j ) + b[l * _nTotal + j];
        }
    }
}

void Closure::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < _nTotal; ++j ) {
            b[l * _nTotal + j] = b[l * _nTotal + j] - A( l, j );
        }
    }
}

void Closure::DS( Vector& ds, const Vector& u ) const {}

std::vector<Polynomial*> Closure::GetBasis() { return _basis; }

std::vector<Polynomial*> Closure::GetQuadrature() { return _quad; }

void Closure::SetAlpha( double alpha ) { _alpha = alpha; }

void Closure::SolveClosureSafe( Matrix& lambda, const Matrix& u ) {
    int maxRefinements = 1000;

    Matrix H( _nStates * _nTotal, _nStates * _nTotal );
    Vector g( _nStates * _nTotal );
    Vector dlambdaNew( _nStates * _nTotal );

    Vector dlambda = -g;
    Matrix lambdaNew( _nStates, _nTotal );

    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        Gradient( g, lambda, u );
        dlambda = -g;
        Hessian( H, lambda );
        posv( H, g );
        AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew );
        Gradient( dlambdaNew, lambdaNew, u );
        int refinementCounter = 0;
        while( CalcNorm( dlambda ) < CalcNorm( dlambdaNew ) || !std::isfinite( CalcNorm( dlambdaNew ) ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew );
            Gradient( dlambdaNew, lambdaNew, u );
            if( CalcNorm( dlambdaNew ) < _settings->GetEpsilon() ) {
                lambda = lambdaNew;
                return;
            }
            else if( ++refinementCounter > maxRefinements ) {
                _log->error( "[closure] Newton needed too many refinement steps!" );
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( CalcNorm( dlambdaNew ) < _settings->GetEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

void Closure::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }

unsigned Closure::GetMaxIterations() const { return _maxIterations; }
