#include "closure.h"
#include "boundedbarrier.h"
#include "eulerclosure.h"
#include "eulerclosure2d.h"
#include "stochasticgalerkin.h"

Closure::Closure( Settings* settings )
    : _settings( settings ), _nMoments( _settings->GetNMoments() ), _nQuadPoints( _settings->GetNQuadPoints() ), _nStates( _settings->GetNStates() ) {
    _log = spdlog::get( "event" );
    // initialize classes
    _basis = Polynomial::Create( _settings, _nMoments );
    _quad  = Polynomial::Create( _settings, _nQuadPoints );

    // compute ntotal number of moments and quad points
    _numDimXi = _settings->GetNDimXi();
    _nQTotal  = std::pow( _settings->GetNQuadPoints(), _numDimXi );
    _nTotal   = std::pow( _settings->GetNMoments(), _numDimXi );

    // calculate basis functions evaluated at the quadrature points
    _phiTilde    = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeWf  = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeVec = std::vector<Vector>( _nQTotal, Vector( _nTotal, 0.0 ) );

    // setup map from k total to individual indices
    std::vector<std::vector<unsigned>> indices;
    std::vector<std::vector<unsigned>> indicesQ;
    indicesQ.resize( _nQTotal );
    indices.resize( _nTotal );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        indicesQ[k].resize( _numDimXi );
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indicesQ[k][l] = unsigned( ( k - k % unsigned( std::pow( _nQuadPoints, l ) ) ) / unsigned( std::pow( _nQuadPoints, l ) ) ) % _nQuadPoints;
            // std::cout << indicesQ[k][l] << " ";
        }
        // std::cout << std::endl;
    }

    for( unsigned i = 0; i < _nTotal; ++i ) {
        indices[i].resize( _numDimXi );
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indices[i][l] = unsigned( ( i - i % unsigned( std::pow( _nMoments, l ) ) ) / unsigned( std::pow( _nMoments, l ) ) ) % _nMoments;
            // std::cout << indices[i][l] << " ";
        }
        // std::cout << std::endl;
    }

    Vector xi = _quad->GetNodes();
    Vector w  = _quad->GetWeights();

    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned l = 0; l < _numDimXi; ++l ) {
                // for( unsigned n = 0; n < _numDimXi; ++n ) {
                _phiTilde( k, i ) *=
                    _basis->Evaluate( indices[i][l], xi[indicesQ[k][l]] ) / _basis->L2Norm( indices[i][l] );    // sqrt( 2.0 * i + 1.0 );
                _phiTildeWf( k, i ) *= _basis->Evaluate( indices[i][l], xi[indicesQ[k][l]] ) / _basis->L2Norm( indices[i][l] ) * w[indicesQ[k][l]] *
                                       _basis->fXi( xi[indicesQ[k][l]] );
                //}
            }
            // multiplied by pdf
            _phiTildeVec[k][i] = _phiTilde( k, i );    // sqrt( 2.0 * i + 1.0 );
        }
    }

    // Test

    // std::cout << _phiTilde << std::endl;
    // exit( EXIT_FAILURE );
    _phiTildeTrans       = trans( _phiTilde );
    auto phiTildeWfTrans = trans( _phiTildeWf );
    // calculate partial matrix for Hessian calculation
    _hPartial = MatVec( _nQTotal, Matrix( _nTotal, _nTotal ) );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        _hPartial[k] = outer( column( _phiTildeTrans, k ), column( phiTildeWfTrans, k ) );    // TODO
    }
}

Closure::~Closure() {
    delete _basis;
    delete _quad;
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
    else {
        log->error( "Invalid closure type" );
        exit( EXIT_FAILURE );
    }
}

void Closure::SolveClosure( Matrix& lambda, const Matrix& u ) {
    // std::cout << "SolveClosure started..." << std::endl;
    int maxRefinements = 1000;

    Matrix H( _nStates * _nTotal, _nStates * _nTotal, 0.0 );
    Vector g( _nStates * _nTotal, 0.0 );
    Vector dlambdaNew( _nStates * _nTotal, 0.0 );

    // check if initial guess is good enough
    Gradient( g, lambda, u );
    if( CalcNorm( g ) < _settings->GetEpsilon() ) {
        return;
    }
    // std::cout << "before first Hessian inversion..." << std::endl;
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda );
    posv( H, g );
    // std::cout << "... done!" << std::endl;
    if( _settings->GetMaxIterations() == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda );
        return;
    }
    Matrix lambdaNew( _nStates, _nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew );
    Gradient( dlambdaNew, lambdaNew, u );
    // perform Newton iterations
    for( unsigned l = 0; l < _settings->GetMaxIterations(); ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, u );
            dlambda = -g;
            Hessian( H, lambda );
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew );
            Gradient( dlambdaNew, lambdaNew, u );
        }
        int refinementCounter = 0;
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

void Closure::EvaluateLambda( Matrix& out, const Matrix& lambda ) const { out = lambda * _phiTildeTrans; }

Vector Closure::EvaluateLambda( const Matrix& lambda, const Vector& xi, unsigned k ) {
    Vector out = Vector( _nStates, 0.0 );
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            out[l] += lambda( l, i ) * _basis->Evaluate( i, xi[k] ) * ( 2.0 * i + 1.0 );
        }
    }
    return out;
}

// TODO
Matrix Closure::EvaluateLambda( const Matrix& lambda, const Vector& xi ) {
    Matrix out( _nStates, xi.size(), 0.0 );
    for( unsigned k = 0; k < xi.size(); ++k ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned i = 0; i < _nTotal; ++i ) {
                out( l, k ) += lambda( l, i ) * _basis->Evaluate( i, xi[k] ) * ( 2.0 * i + 1.0 );
            }
        }
    }
    return out;
}

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u ) {
    Vector uKinetic( _nStates, 0.0 );
    g.reset();

    for( unsigned k = 0; k < _nQTotal; ++k ) {
        U( uKinetic, lambda * _phiTildeVec[k] );
        for( unsigned j = 0; j < _nTotal; ++j ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * _nTotal + j] += uKinetic[l] * _phiTildeWf( k, j );
            }
        }
    }

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
