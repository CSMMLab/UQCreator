#include "thermalradiationclosure.h"

ThermalRadiationClosure::ThermalRadiationClosure( Settings* settings ) : Closure( settings ), _nHydroStates( 1 ), _nMoments( 2 ) { _alpha = 1.0; }

ThermalRadiationClosure::~ThermalRadiationClosure() {}

void ThermalRadiationClosure::U( Vector& out, const Vector& Lambda ) {
    out[0] = Lambda[0];
    out[1] = Lambda[1];
    out[2] = exp( Lambda[2] );
}

void ThermalRadiationClosure::U( Vector& out, const Vector& Lambda, bool dummy ) { out[0] = exp( Lambda[0] ); }

void ThermalRadiationClosure::U( Tensor& out, const Tensor& Lambda ) {
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            out( 0, l, k ) = Lambda( 0, l, k );
            out( 1, l, k ) = Lambda( 1, l, k );
            out( 2, l, k ) = exp( Lambda( 2, l, k ) );
        }
    }
}

Tensor ThermalRadiationClosure::U( const Tensor& Lambda ) {
    Tensor y( _nStates, _nMultiElements, Lambda.columns(), 0.0 );
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        for( unsigned k = 0; k < Lambda.columns(); ++k ) {
            y( 0, l, k ) = Lambda( 0, l, k );
            y( 1, l, k ) = Lambda( 1, l, k );
            y( 2, l, k ) = exp( Lambda( 2, l, k ) );
        }
    }

    return y;
}

void ThermalRadiationClosure::DU( Matrix& y, const Vector& Lambda ) { y( 0, 0 ) = exp( Lambda[0] ); }

void ThermalRadiationClosure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector uKinetic( _nHydroStates, 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    g.reset();

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ), true );
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nHydroStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeF( k, i ) * _wGrid[refLevel][k];
            }
        }
    }

    SubstractVectorMatrixOnVector( g, u, _nTotalForRef[refLevel] );
}

void ThermalRadiationClosure::Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel ) {
    H.reset();
    Matrix dUdLambda( _nHydroStates, _nHydroStates );    // TODO: preallocate Matrix for Hessian computation -> problems omp
    unsigned nTotal = _nTotalForRef[refLevel];

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dUdLambda, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned l = 0; l < _nHydroStates; ++l ) {
            for( unsigned m = 0; m < _nHydroStates; ++m ) {
                for( unsigned j = 0; j < nTotal; ++j ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        H( m * nTotal + j, l * nTotal + i ) += _hPartial[k]( j, i ) * _wGrid[refLevel][k] * dUdLambda( l, m );
                    }
                }
            }
        }
    }
}

void ThermalRadiationClosure::SolveClosure( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel ) {
    int maxRefinements = 1000;
    unsigned nTotal    = _nTotalForRef[refLevel];

    Vector g( _nHydroStates * nTotal );

    // save hydro part on lambda and u
    Matrix lambda( _nHydroStates, nTotal );
    Matrix u( _nHydroStates, nTotal );
    for( unsigned s = 0; s < _nHydroStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            lambda( s, i ) = lambdaFull( _nMoments + s, i );
            u( s, i )      = uFull( _nMoments + s, i );
        }
    }
    // std::cout << u << std::endl;
    // std::cout << lambda << std::endl;

    // save SG result for radiation part
    for( unsigned s = 0; s < _nMoments; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            lambdaFull( s, i ) = uFull( s, i );
        }
    }

    // check if initial guess is good enough
    Gradient( g, lambda, u, refLevel );
    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
        return;
    }
    Matrix H( _nHydroStates * nTotal, _nHydroStates * nTotal );
    Vector dlambdaNew( _nHydroStates * nTotal );
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda, refLevel );
    // std::cout << "H = " << std::endl << H << std::endl;
    // std::cout << "g = " << std::endl << g << std::endl;
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
        for( unsigned s = 0; s < _nHydroStates; ++s ) {
            for( unsigned i = 0; i < nTotal; ++i ) {
                lambdaFull( _nMoments + s, i ) = lambda( s, i );
            }
        }
        return;
    }
    Matrix lambdaNew( _nHydroStates, nTotal );
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
        // std::cout << "H = " << std::endl << H << std::endl;
        // std::cout << "g = " << std::endl << g << std::endl;
        // std::cout << "Res " << CalcNorm( dlambdaNew, nTotal ) << std::endl;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, u, refLevel );
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                for( unsigned s = 0; s < _nHydroStates; ++s ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        lambdaFull( _nMoments + s, i ) = lambdaNew( s, i );
                    }
                }
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
            for( unsigned s = 0; s < _nHydroStates; ++s ) {
                for( unsigned i = 0; i < nTotal; ++i ) {
                    lambdaFull( _nMoments + s, i ) = lambdaNew( s, i );
                }
            }
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

void ThermalRadiationClosure::SolveClosureSafe( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel ) {
    this->SolveClosure( lambdaFull, uFull, refLevel );
}

void ThermalRadiationClosure::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            y( l, j ) = A( l, j ) + b[l * nTotal + j];
        }
    }
}

void ThermalRadiationClosure::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            b[l * nTotal + j] = b[l * nTotal + j] - A( l, j );
        }
    }
}

double ThermalRadiationClosure::CalcNorm( Vector& test, unsigned nTotal ) const {
    double out = 0.0;
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out += pow( test[l * nTotal + i], 2 );
        }
    }
    return sqrt( out );
}
