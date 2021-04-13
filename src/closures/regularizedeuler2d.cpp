#include "closures/regularizedeuler2d.h"

RegularizedEuler2D::RegularizedEuler2D( Settings* settings ) : EulerClosure2D( settings ), _eta( _settings->GetRegularizationStrength() ) {}

RegularizedEuler2D::~RegularizedEuler2D() {}

void RegularizedEuler2D::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector uKinetic( _nStates, 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    g.reset();

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeF( k, i ) * _wGrid[refLevel][k];
            }
        }
    }

    for( unsigned i = 0; i < nTotal; ++i ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            g[l * nTotal + i] += _eta * lambda( l, i );
        }
    }

    SubstractVectorMatrixOnVector( g, u, _nTotalForRef[refLevel] );
}

void RegularizedEuler2D::GradientNoRegularizaton( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector uKinetic( _nStates, 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    g.reset();

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeF( k, i ) * _wGrid[refLevel][k];
            }
        }
    }

    SubstractVectorMatrixOnVector( g, u, _nTotalForRef[refLevel] );
}

void RegularizedEuler2D::Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel ) {
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
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            H( l * nTotal + i, l * nTotal + i ) += _eta;
        }
    }
}

void RegularizedEuler2D::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
    int maxRefinements = 1000;
    unsigned nTotal    = _nTotalForRef[refLevel];
    Matrix lambdaMat   = Matrix( _settings->GetNStates(), nTotal );
    Matrix uMat        = Matrix( _settings->GetNStates(), nTotal );
    bool breakFlag     = false;

    Vector g( _nStates * nTotal );
    Matrix H( _nStates * nTotal, _nStates * nTotal );
    for( unsigned l = 0; l < _nMultiElements; ++l ) {
        // save solution in element l as a matrix
        _filter->FilterMoments( uMat, u, l );
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned i = 0; i < nTotal; ++i ) {
                lambdaMat( s, i ) = lambda( s, l, i );
            }
        }

        // check if initial guess is good enough
        // check if initial guess is good enough
        GradientNoRegularizaton( g, lambdaMat, uMat, refLevel );
        if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
            return;
        }

        Gradient( g, lambdaMat, uMat, refLevel );

        Vector dlambdaNew( _nStates * nTotal );
        // calculate initial Hessian and gradient
        Vector dlambda = -g;
        Hessian( H, lambdaMat, refLevel );

        posv( H, g );
        if( _maxIterations == 1 ) {
            AddMatrixVectorToMatrix( lambdaMat, -_alpha * g, lambdaMat, nTotal );
            for( unsigned s = 0; s < _nStates; ++s ) {
                for( unsigned i = 0; i < nTotal; ++i ) {
                    lambda( s, l, i ) = lambdaMat( s, i );
                }
            }
            continue;
        }
        Matrix lambdaNew( _nStates, nTotal );
        AddMatrixVectorToMatrix( lambdaMat, -_alpha * g, lambdaNew, nTotal );
        Gradient( dlambdaNew, lambdaNew, uMat, refLevel );
        // perform Newton iterations
        unsigned m;
        for( m = 0; m < _maxIterations; ++m ) {
            double stepSize = 1.0;
            if( m != 0 ) {
                Gradient( g, lambdaMat, uMat, refLevel );
                dlambda = -g;
                Hessian( H, lambdaMat, refLevel );
                posv( H, g );
                AddMatrixVectorToMatrix( lambdaMat, -stepSize * _alpha * g, lambdaNew, nTotal );
                Gradient( dlambdaNew, lambdaNew, uMat, refLevel );
            }
            int refinementCounter = 0;
            // std::cout << CalcNorm( dlambdaNew, nTotal ) << std::endl;
            while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
                stepSize *= 0.5;
                AddMatrixVectorToMatrix( lambdaMat, -stepSize * _alpha * g, lambdaNew, nTotal );
                Gradient( dlambdaNew, lambdaNew, uMat, refLevel );
                if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                    for( unsigned s = 0; s < _nStates; ++s ) {
                        for( unsigned i = 0; i < nTotal; ++i ) {
                            lambda( s, l, i ) = lambdaNew( s, i );
                        }
                    }
                    breakFlag = true;
                    break;
                }
                else if( ++refinementCounter > maxRefinements ) {
                    _log->error( "[closure] Newton needed too many refinement steps!" );
                    exit( EXIT_FAILURE );
                }
            }
            if( breakFlag ) {
                breakFlag = false;
                break;
            }
            lambdaMat = lambdaNew;
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                for( unsigned s = 0; s < _nStates; ++s ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        lambda( s, l, i ) = lambdaNew( s, i );
                    }
                }
                break;
            }
        }
        if( m == _maxIterations ) {
            _log->error( "[closure] Newton did not converge!" );
            exit( EXIT_FAILURE );
        }
    }
}

void RegularizedEuler2D::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) { SolveClosure( lambda, u, refLevel ); }
