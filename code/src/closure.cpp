#include "closure.h"
#include "boundedbarrier.h"
#include "stochasticgalerkin.h"

Closure::Closure( Problem* problem )
    : _problem( problem ), _nMoments( _problem->GetNMoments() ), _nQuadPoints( _problem->GetNQuadPoints() ), _nStates( _problem->GetNStates() ) {
    // initialize classes
    _basis = new Legendre( _nMoments );
    _quad  = new Legendre( _nQuadPoints );
    // calculate basis functions evaluated at the quadrature points
    _phi         = std::vector<blaze::DynamicVector<double>>( _nQuadPoints, blaze::DynamicVector<double>( _nMoments, 0.0 ) );
    _phiTilde    = blaze::DynamicMatrix<double>( _nQuadPoints, _nMoments, 0.0 );
    _phiTildeW   = blaze::DynamicMatrix<double>( _nQuadPoints, _nMoments, 0.0 );
    _phiTildeVec = std::vector<blaze::DynamicVector<double>>( _nQuadPoints, blaze::DynamicVector<double>( _nMoments, 0.0 ) );

    blaze::DynamicVector<double> xi = _quad->GetNodes();
    blaze::DynamicVector<double> w  = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ) {
        for( int i = 0; i < _nMoments; ++i ) {
            _phi[k][i]         = _basis->Evaluate( i, xi[k] );
            _phiTilde( k, i )  = _phi[k][i] * ( 2.0 * i + 1.0 );
            _phiTildeW( k, i ) = _phiTilde( k, i ) * w[k];
            _phiTildeVec[k][i] = _phi[k][i] * ( 2.0 * i + 1.0 );
        }
    }
    _phiTildeTrans = blaze::trans( _phiTilde );
    // calculate partial matrix for Hessian calculation
    _hPartial = std::vector<blaze::DynamicMatrix<double>>( _nQuadPoints, blaze::DynamicMatrix<double>( _nMoments, _nMoments, 0.0 ) );
    for( int k = 0; k < _nQuadPoints; ++k ) {
        _hPartial[k] = blaze::trans( row( _phiTilde, k ) ) * row( _phiTilde, k ) * w[k];
    }
    double du = 0.0;
    _uMinus   = 3.0 - du;
    _uPlus    = 12.0 + du;
}

Closure* Closure::Create( Problem* problem ) {
    std::string closure = problem->GetClosureType();
    if( closure.compare( "BoundedBarrier" ) == 0 || closure.compare( "BB" ) == 0 ) {
        return new BoundedBarrier( problem );
    }
    else if( closure.compare( "StochasticGalerkin" ) == 0 || closure.compare( "SG" ) == 0 ) {
        return new StochasticGalerkin( problem );
    }
    else {
        std::cerr << "Invalid closure type" << std::endl;
        exit( EXIT_FAILURE );
        return NULL;
    }
}

blaze::DynamicMatrix<double> Closure::SolveClosure( const blaze::DynamicMatrix<double>& u, blaze::DynamicMatrix<double>& lambda ) {
    int maxRefinements = 1000;
    blaze::DynamicMatrix<double> H( _nStates * _nMoments, _nStates * _nMoments, 0.0 );
    blaze::DynamicVector<double> g( _nStates * _nMoments, 0.0 );
    blaze::DynamicVector<double> dlambdaNew( _nStates * _nMoments, 0.0 );

    // check if initial guess is good enough
    Gradient( g, lambda, u );

    if( CalcNorm( g ) < _problem->GetEpsilon() ) {
        return lambda;
    }
    // calculate initial Hessian and gradient
    blaze::DynamicVector<double> dlambda = -g;
    Hessian( H, lambda );
    blaze::posv( H, g, 'L' );
    blaze::DynamicMatrix<double> lambdaNew = lambda - MakeMatrix( g );
    Gradient( dlambdaNew, lambdaNew, u );
    // perform Newton iterations
    for( int l = 0; l < _problem->GetMaxIterations(); ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            Gradient( g, lambda, u );
            dlambda = -g;
            Hessian( H, lambda );
            blaze::posv( H, g, 'L' );
            lambdaNew = lambda - stepSize * MakeMatrix( g );
            Gradient( dlambdaNew, lambdaNew, u );
        }
        int refinementCounter = 0;
        while( CalcNorm( dlambda ) < CalcNorm( dlambdaNew ) ) {
            stepSize *= 0.5;
            lambdaNew = lambda - stepSize * MakeMatrix( g );
            Gradient( dlambdaNew, lambdaNew, u );
            if( CalcNorm( dlambdaNew ) < _problem->GetEpsilon() ) {
                return lambdaNew;
            }
            else if( ++refinementCounter > maxRefinements ) {
                std::cerr << "[ERROR]: Newton needed too many refinement steps!" << std::endl;
                exit( EXIT_FAILURE );
            }
        }
        lambda = lambdaNew;

        if( CalcNorm( dlambdaNew ) < _problem->GetEpsilon() ) {
            return lambdaNew;
        }
    }
    std::cerr << "[ERROR]: Newton did not converge!" << std::endl;
    exit( EXIT_FAILURE );
    return blaze::DynamicMatrix<double>( _nStates, _nMoments, -1.0 );
}

double Closure::CalcNorm( blaze::DynamicVector<double>& test ) {
    double out = 0.0;
    double tmp;
    for( int l = 0; l < _nStates; ++l ) {
        tmp = 0.0;
        for( int i = 0; i < _nMoments; ++i ) {
            tmp += pow( test[l * _nMoments + i], 2 );
        }
        out += sqrt( tmp );
    }
    return out;
}

blaze::DynamicVector<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, int k ) { return lambda * _phiTildeVec[k]; }

blaze::DynamicMatrix<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda ) const { return lambda * _phiTildeTrans; }

blaze::DynamicVector<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicVector<double>& xi, int k ) {
    blaze::DynamicVector<double> out = blaze::DynamicVector<double>( _nStates, 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        for( int i = 0; i < _nMoments; ++i ) {
            out[l] += lambda( l, i ) * _basis->Evaluate( i, xi[k] ) * ( 2.0 * i + 1.0 );
        }
    }
    return out;
}

// TODO
blaze::DynamicMatrix<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicVector<double>& xi ) {
    blaze::DynamicMatrix<double> out( _nStates, xi.size(), 0.0 );
    for( unsigned int k = 0; k < xi.size(); ++k ) {
        for( int l = 0; l < _nStates; ++l ) {
            for( int i = 0; i < _nMoments; ++i ) {
                out( l, k ) += lambda( l, i ) * _basis->Evaluate( i, xi[k] ) * ( 2.0 * i + 1.0 );
            }
        }
    }
    return out;
}

void Closure::Gradient( blaze::DynamicVector<double>& g, const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicMatrix<double>& u ) {
    blaze::DynamicVector<double> uKinetic( _nStates, 0.0 );
    g.reset();

    for( int k = 0; k < _nQuadPoints; ++k ) {
        U( uKinetic, lambda * _phiTildeVec[k] );
        for( int j = 0; j < _nMoments; ++j ) {
            for( int l = 0; l < _nStates; ++l ) {
                g[l * _nMoments + j] += uKinetic[l] * _phiTildeW( k, j );
            }
        }
    }

    g = 0.5 * g - MakeVector( u );
}

void Closure::Hessian( blaze::DynamicMatrix<double>& H, const blaze::DynamicMatrix<double>& lambda ) {
    blaze::DynamicMatrix<double> dLambdadU( _nStates, _nStates, 0.0 );
    H.reset();

    for( int k = 0; k < _nQuadPoints; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dLambdadU, lambda * _phiTildeVec[k] );
        for( int l = 0; l < _nStates; ++l ) {
            for( int m = 0; m < _nStates; ++m ) {
                for( int j = 0; j < _nMoments; ++j ) {
                    for( int i = 0; i < _nMoments; ++i ) {
                        H( m * _nMoments + j, l * _nMoments + i ) += 0.5 * _hPartial[k]( j, i ) * dLambdadU( l, m );
                    }
                }
            }
        }
    }
}

blaze::DynamicVector<double> Closure::MakeVector( const blaze::DynamicMatrix<double>& mat ) const {
    blaze::DynamicVector<double> y( _nStates * _nMoments, 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        for( int j = 0; j < _nMoments; ++j ) {
            y[l * _nMoments + j] = mat( l, j );
        }
    }
    return y;
}

blaze::DynamicMatrix<double> Closure::MakeMatrix( const blaze::DynamicVector<double>& vec ) const {
    blaze::DynamicMatrix<double> y( _nStates, _nMoments, 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        for( int j = 0; j < _nMoments; ++j ) {
            y( l, j ) = vec[l * _nMoments + j];
        }
    }
    return y;
}
