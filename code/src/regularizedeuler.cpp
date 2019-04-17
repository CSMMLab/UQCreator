#include "regularizedeuler.h"

RegularizedEuler::RegularizedEuler( Settings* settings ) : EulerClosure2D( settings ), _eta( 1e-9 ), _lambda( 0.00015 ) {
    unsigned nMoments = _settings->GetNMoments();
    _filterFunction   = Vector( _settings->GetNTotal(), 1.0 );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                // if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                // if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                unsigned index = unsigned( ( i - i % unsigned( std::pow( nMoments, l ) ) ) / unsigned( std::pow( nMoments, l ) ) ) % nMoments;
                _filterFunction[i] *= 1.0 / ( 1.0 + _lambda * pow( index, 2 ) * pow( index + 1, 2 ) );
            }
        }
    }
}

RegularizedEuler::~RegularizedEuler() {}

void RegularizedEuler::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal ) {
    Vector uKinetic( _nStates, 0.0 );
    g.reset();

    for( unsigned k = 0; k < nQTotal; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned l = 0; l < _nStates; ++l ) {
                g[l * nTotal + i] += uKinetic[l] * _phiTildeWf( k, i );
            }
        }
    }

    for( unsigned i = 0; i < nTotal; ++i ) {
        for( unsigned l = 0; l < _nStates; ++l ) {
            g[l * nTotal + i] += _eta * lambda( l, i );
        }
    }

    SubstractVectorMatrixOnVector( g, u, nTotal );
}

void RegularizedEuler::Hessian( Matrix& H, const Matrix& lambda, unsigned nTotal, unsigned nQTotal ) {
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
    for( unsigned l = 0; l < _nStates; ++l ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            H( l * nTotal + i, l * nTotal + i ) += _eta;
        }
    }
}

void RegularizedEuler::SolveClosure( Matrix& lambda, const Matrix& u, unsigned nTotal, unsigned nQTotal ) {
    Matrix uF( _settings->GetNStates(), nTotal );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            uF( s, i ) = _filterFunction[i] * u( s, i );
        }
    }
    int maxRefinements = 1000;

    Vector g( _nStates * nTotal );

    // check if initial guess is good enough
    Gradient( g, lambda, uF, nTotal, nQTotal );
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
    Gradient( dlambdaNew, lambdaNew, uF, nTotal, nQTotal );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, uF, nTotal, nQTotal );
            dlambda = -g;
            Hessian( H, lambda, nTotal, nQTotal );
            // std::cout << H << std::endl;
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, nTotal, nQTotal );
        }
        int refinementCounter = 0;
        // std::cout << "Res is " << CalcNorm( dlambda, nTotal ) << std::endl;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, nTotal, nQTotal );
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
