#include "euler1dfpfilter.h"

Euler1DFPFilter::Euler1DFPFilter( Settings* settings ) : EulerClosure( settings ), _lambda( _settings->GetFilterStrength() ) {
    unsigned maxDegree = _settings->GetMaxDegree();
    _filterFunction    = Vector( _settings->GetNTotal(), 1.0 );
    _hyperbolicityFix  = false;    // turns on Graham's hyperbolicity fix
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < _settings->GetNTotal(); ++i ) {
            for( unsigned l = 0; l < _settings->GetNDimXi(); ++l ) {
                unsigned index =
                    unsigned( ( i - i % unsigned( std::pow( maxDegree + 1, l ) ) ) / unsigned( std::pow( maxDegree + 1, l ) ) ) % ( maxDegree + 1 );
                _filterFunction[i] *= std::exp( -_lambda * index * ( index + 1 ) );
            }
        }
    }
}

Matrix Euler1DFPFilter::Filter( const Matrix& u, unsigned nTotal ) const {
    Matrix uF( _settings->GetNStates(), nTotal );
    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            uF( s, i ) = _filterFunction[i] * u( s, i );
        }
    }
    return uF;
}

void Euler1DFPFilter::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    unsigned nTotal    = _nTotalForRef[refLevel];
    int maxRefinements = 1000;
    // check if initial guess is good enough
    Vector g( _nStates * nTotal );
    Gradient( g, lambda, u, refLevel );
    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
        return;
    }

    Matrix uF = Filter( u, nTotal );

    Gradient( g, lambda, uF, refLevel );

    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda, refLevel );
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
        if( _hyperbolicityFix ) lambda = Filter( lambda, nTotal );
        return;
    }
    Matrix lambdaNew( _nStates, nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew, nTotal );
    Gradient( dlambdaNew, lambdaNew, uF, refLevel );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, uF, refLevel );
            dlambda = -g;
            Hessian( H, lambda, refLevel );
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, refLevel );
        }
        int refinementCounter = 0;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, refLevel );
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                if( _hyperbolicityFix )
                    lambda = Filter( lambdaNew, nTotal );
                else
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
            if( _hyperbolicityFix )
                lambda = Filter( lambdaNew, nTotal );
            else
                lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}

void Euler1DFPFilter::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    unsigned nTotal    = _nTotalForRef[refLevel];
    int maxRefinements = 1000;

    // check if initial guess is good enough
    Vector g( _nStates * nTotal );
    Gradient( g, lambda, u, refLevel );
    if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
        return;
    }

    Matrix uF = Filter( u, nTotal );

    Gradient( g, lambda, uF, refLevel );
    Matrix H( _nStates * nTotal, _nStates * nTotal );
    Vector dlambdaNew( _nStates * nTotal );
    // calculate initial Hessian and gradient
    Vector dlambda = -g;
    // std::cout << g << std::endl;
    Hessian( H, lambda, refLevel );
    posv( H, g );
    if( _maxIterations == 1 ) {
        AddMatrixVectorToMatrix( lambda, -_alpha * g, lambda, nTotal );
        if( _hyperbolicityFix ) lambda = Filter( lambda, nTotal );
        return;
    }
    Matrix lambdaNew( _nStates, nTotal );
    AddMatrixVectorToMatrix( lambda, -_alpha * g, lambdaNew, nTotal );
    Gradient( dlambdaNew, lambdaNew, uF, refLevel );
    // perform Newton iterations
    for( unsigned l = 0; l < _maxIterations; ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, uF, refLevel );
            dlambda = -g;
            Hessian( H, lambda, refLevel );
            posv( H, g );
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, refLevel );
        }
        int refinementCounter = 0;
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) || !std::isfinite( CalcNorm( dlambdaNew, nTotal ) ) ) {
            stepSize *= 0.5;
            AddMatrixVectorToMatrix( lambda, -stepSize * _alpha * g, lambdaNew, nTotal );
            Gradient( dlambdaNew, lambdaNew, uF, refLevel );
            if( CalcNorm( dlambdaNew, nTotal ) < _settings->GetEpsilon() ) {
                if( _hyperbolicityFix )
                    lambda = Filter( lambdaNew, nTotal );
                else
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
            if( _hyperbolicityFix )
                lambda = Filter( lambdaNew, nTotal );
            else
                lambda = lambdaNew;
            return;
        }
    }
    _log->error( "[closure] Newton did not converge!" );
    exit( EXIT_FAILURE );
}
