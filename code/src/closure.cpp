#include "closure.h"
#include "boundedbarrier.h"
#include "eulerclosure.h"
#include "eulerclosure2d.h"
#include "stochasticgalerkin.h"

Closure::Closure( Settings* settings )
    : _settings( settings ), _nMoments( _settings->GetNMoments() ), _nQuadPoints( _settings->GetNQuadPoints() ), _nStates( _settings->GetNStates() ) {
    // initialize classes
    _basis = new Legendre( _nMoments );
    _quad  = new Legendre( _nQuadPoints );
    // calculate basis functions evaluated at the quadrature points
    _phi         = std::vector<Vector>( _nQuadPoints, Vector( _nMoments, 0.0 ) );
    _phiTilde    = Matrix( _nQuadPoints, _nMoments, 0.0 );
    _phiTildeW   = Matrix( _nQuadPoints, _nMoments, 0.0 );
    _phiTildeVec = std::vector<Vector>( _nQuadPoints, Vector( _nMoments, 0.0 ) );

    Vector xi = _quad->GetNodes();
    Vector w  = _quad->GetWeights();

    for( unsigned k = 0; k < _nQuadPoints; ++k ) {
        for( unsigned i = 0; i < _nMoments; ++i ) {
            _phi[k][i]         = _basis->Evaluate( i, xi[k] );
            _phiTilde( k, i )  = _phi[k][i] * ( 2.0 * i + 1.0 );    // sqrt( 2.0 * i + 1.0 );
            _phiTildeW( k, i ) = _phiTilde( k, i ) * w[k];
            _phiTildeVec[k][i] = _phi[k][i] * ( 2.0 * i + 1.0 );    // sqrt( 2.0 * i + 1.0 );
        }
    }
    _phiTildeTrans = trans( _phiTilde );
    // calculate partial matrix for Hessian calculation
    _hPartial = MatVec( _nQuadPoints, Matrix( _nMoments, _nMoments, 0.0 ) );
    for( unsigned k = 0; k < _nQuadPoints; ++k ) {
        _hPartial[k] = outer( column( _phiTildeTrans, k ), column( _phiTildeTrans, k ) ) * w[k];
    }

    _perm = new int[_nStates * _nMoments];
    for( int i = 0; i < static_cast<int>( _nStates * _nMoments ); ++i ) {
        _perm[i] = i;
    }
}

Closure::~Closure() {
    delete _basis;
    delete _quad;
    delete _perm;
}

Closure* Closure::Create( Settings* settings ) {
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
        std::cerr << "Invalid closure type" << std::endl;
        exit( EXIT_FAILURE );
    }
}

void Closure::SolveClosure( Matrix& lambda, const Matrix& u ) {
    int maxRefinements = 1000;

    Matrix H( _nStates * _nMoments, _nStates * _nMoments, 0.0 );
    Vector g( _nStates * _nMoments, 0.0 );
    Vector dlambdaNew( _nStates * _nMoments, 0.0 );

    // check if initial guess is good enough
    Gradient( g, lambda, u );
    if( CalcNorm( g ) < _settings->GetEpsilon() ) {
        return;
    }
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    Hessian( H, lambda );
    posv( H, g );
    // gesv( H, g, _perm );
    Matrix lambdaNew( _nStates, _nMoments );
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
            // gesv( H, g, _perm );
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
                std::cerr << "[ERROR]: Newton needed too many refinement steps!" << std::endl;
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;
        if( CalcNorm( dlambdaNew ) < _settings->GetEpsilon() ) {
            lambda = lambdaNew;
            return;
        }
    }
    std::cerr << "[Closure][ERROR]: Newton did not converge!" << std::endl;
    exit( EXIT_FAILURE );
}

double Closure::CalcNorm( Vector& test ) {
    double out = 0.0;
    double tmp;
    for( unsigned l = 0; l < _nStates; ++l ) {
        tmp = 0.0;
        for( unsigned i = 0; i < _nMoments; ++i ) {
            tmp += pow( test[l * _nMoments + i], 2 );
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
        for( unsigned i = 0; i < _nMoments; ++i ) {
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
            for( unsigned i = 0; i < _nMoments; ++i ) {
                out( l, k ) += lambda( l, i ) * _basis->Evaluate( i, xi[k] ) * ( 2.0 * i + 1.0 );
            }
        }
    }
    return out;
}

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u ) {
    Vector uKinetic( _nStates, 0.0 );
    g.reset();

    for( unsigned k = 0; k < _nQuadPoints; ++k ) {
        U( uKinetic, lambda * _phiTildeVec[k] );
        for( unsigned j = 0; j < _nMoments; ++j ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * _nMoments + j] += 0.5 * uKinetic[l] * _phiTildeW( k, j );
            }
        }
    }

    SubstractVectorMatrixOnVector( g, u );
}

void Closure::Hessian( Matrix& H, const Matrix& lambda ) {
    Matrix dLambdadU( _nStates, _nStates, 0.0 );
    H.reset();

    for( unsigned k = 0; k < _nQuadPoints; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dLambdadU, lambda * _phiTildeVec[k] );
        for( unsigned l = 0; l < _nStates; ++l ) {
            for( unsigned m = 0; m < _nStates; ++m ) {
                for( unsigned j = 0; j < _nMoments; ++j ) {
                    for( unsigned i = 0; i < _nMoments; ++i ) {
                        H( m * _nMoments + j, l * _nMoments + i ) += 0.5 * _hPartial[k]( j, i ) * dLambdadU( l, m );
                    }
                }
            }
        }
    }
}

void Closure::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < _nMoments; ++j ) {
            y( l, j ) = A( l, j ) + b[l * _nMoments + j];
        }
    }
}

void Closure::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A ) const {
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned j = 0; j < _nMoments; ++j ) {
            b[l * _nMoments + j] = b[l * _nMoments + j] - A( l, j );
        }
    }
}
