#include "closure.h"

Closure::Closure( Problem* problem )
    : _problem( problem ), _nMoments( _problem->GetNMoments() ), _nQuadPoints( _problem->GetNQuadPoints() ), _nStates( _problem->GetNStates() ) {
    // initialize classes
    _basis = new Legendre( _nMoments );
    _quad  = new Legendre( _nQuadPoints );
    // calculate basis functions evaluated at the quadrature points
    _phi                            = std::vector<blaze::DynamicVector<double>>( _nQuadPoints, blaze::DynamicVector<double>( _nMoments, 0.0 ) );
    _phiTilde                       = std::vector<blaze::DynamicVector<double>>( _nQuadPoints, blaze::DynamicVector<double>( _nMoments, 0.0 ) );
    blaze::DynamicVector<double> xi = _quad->GetNodes();
    for( int k = 0; k < _nQuadPoints; ++k ) {
        for( int i = 0; i < _nMoments; ++i ) {
            _phi[k][i]      = _basis->Evaluate( i, xi[k] );
            _phiTilde[k][i] = _phi[k][i] * ( 2.0 * i + 1.0 );
        }
    }
    // calculate partial matrix for Hessian calculation
    _hPartial = std::vector<blaze::DynamicMatrix<double>>( _nQuadPoints, blaze::DynamicMatrix<double>( _nMoments, _nMoments, 0.0 ) );
    blaze::DynamicVector<double> w = _quad->GetWeights();
    for( int k = 0; k < _nQuadPoints; ++k ) {
        for( int i = 0; i < _nMoments; ++i ) {
            for( int j = 0; j < _nMoments; ++j ) {
                _hPartial[k] = _phiTilde[k] * blaze::trans( _phiTilde[k] ) * w[k];
            }
        }
    }
    double du = 0.0;
    _uMinus   = 3.0 - du;
    _uPlus    = 12.0 + du;
}

blaze::DynamicVector<double> Closure::UKinetic( const blaze::DynamicVector<double>& Lambda ) {
    blaze::DynamicVector<double> ePos( _nStates );
    blaze::DynamicVector<double> eNeg( _nStates );
    blaze::DynamicVector<double> out = blaze::DynamicVector<double>( _nStates, 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        ePos[l] = exp( Lambda[l] );
        eNeg[l] = 1 / ePos[l];
        if( Lambda[l] > 0 ) {
            out[l] = _uPlus / ( eNeg[l] + 1.0 ) + _uMinus * eNeg[l] / ( 1.0 + eNeg[l] );
        }
        else {
            out[l] = _uMinus / ( ePos[l] + 1.0 ) + _uPlus * ePos[l] / ( 1.0 + ePos[l] );
        }
    }
    return out;
}

blaze::DynamicMatrix<double> Closure::UKinetic( const blaze::DynamicMatrix<double>& Lambda ) {
    blaze::DynamicMatrix<double> y( _nStates, Lambda.columns(), 0.0 );
    blaze::DynamicMatrix<double> ePos( _nStates, Lambda.columns() );
    blaze::DynamicMatrix<double> eNeg( _nStates, Lambda.columns() );
    for( int l = 0; l < _nStates; ++l ) {
        for( int k = 0; k < Lambda.columns(); ++k ) {
            ePos( l, k ) = exp( Lambda( l, k ) );
            eNeg( l, k ) = 1.0 / ePos( l, k );
            if( Lambda( l, k ) > 0 ) {
                y( l, k ) = _uPlus / ( eNeg( l, k ) + 1.0 ) + _uMinus * eNeg( l, k ) / ( 1.0 + eNeg( l, k ) );
            }
            else {
                y( l, k ) = _uMinus / ( ePos( l, k ) + 1.0 ) + _uPlus * ePos( l, k ) / ( 1.0 + ePos( l, k ) );
            }
        }
    }
    return y;
}

blaze::DynamicMatrix<double> Closure::DUKinetic( const blaze::DynamicVector<double>& Lambda ) {
    blaze::DynamicMatrix<double> y( _nStates, _nStates );
    for( int l = 0; l < _nStates; ++l ) {
        for( int m = 0; m < _nStates; ++m ) {
            double ePos = exp( Lambda[l] );
            double eNeg = 1 / ePos;
            if( Lambda[l] > 0 ) {
                y( l, m ) = eNeg * ( _uPlus - _uMinus ) / ( exp( -2.0 * Lambda[l] ) + 2.0 * eNeg + 1.0 );
            }
            else {
                y( l, m ) = ePos * ( _uPlus - _uMinus ) / ( 1.0 + 2.0 * ePos + exp( 2.0 * Lambda[l] ) );
            }
        }
    }
    return y;
}

blaze::DynamicMatrix<double> Closure::SolveClosure( const blaze::DynamicMatrix<double>& u, blaze::DynamicMatrix<double>& lambda ) {
    int maxRefinements = 1000;

    // check if initial guess is good enough
    blaze::DynamicVector<double> g = Gradient( lambda, u );

    if( CalcNorm( g ) < _problem->GetEpsilon() ) {
        return lambda;
    }
    // calculate initial Hessian and gradient
    blaze::DynamicVector<double> dlambda = -g;
    blaze::DynamicMatrix<double> H       = Hessian( lambda );
    blaze::posv( H, g, 'L' );
    blaze::DynamicMatrix<double> lambdaNew  = lambda - MakeMatrix( g );
    blaze::DynamicVector<double> dlambdaNew = Gradient( lambdaNew, u );
    // perform Newton iterations
    for( int l = 0; l < _problem->GetMaxIterations(); ++l ) {
        double stepSize = 1.0;
        if( l != 0 ) {
            g       = Gradient( lambda, u );
            dlambda = -g;
            H       = Hessian( lambda );
            blaze::posv( H, g, 'L' );
            lambdaNew  = lambda - stepSize * MakeMatrix( g );
            dlambdaNew = Gradient( lambdaNew, u );
        }
        int refinementCounter = 0;
        while( CalcNorm( dlambda ) < CalcNorm( dlambdaNew ) ) {
            stepSize *= 0.5;
            lambdaNew  = lambda - stepSize * MakeMatrix( g );
            dlambdaNew = Gradient( lambdaNew, u );
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

blaze::DynamicVector<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda, int k ) { return lambda * _phiTilde[k]; }

blaze::DynamicMatrix<double> Closure::EvaluateLambda( const blaze::DynamicMatrix<double>& lambda ) const {
    blaze::DynamicMatrix<double> tmp( _nStates, _nQuadPoints, 0.0 );
    for( int l = 0; l < _nStates; ++l ) {
        for( int k = 0; k < _nQuadPoints; ++k ) {
            for( int i = 0; i < _nMoments; ++i ) {
                tmp( l, k ) += lambda( l, i ) * _phiTilde[k][i];
            }
        }
    }
    return tmp;
}

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

blaze::DynamicVector<double> Closure::Gradient( const blaze::DynamicMatrix<double>& lambda, const blaze::DynamicMatrix<double>& u ) {
    blaze::DynamicVector<double> g( _nStates * _nMoments, 0.0 );
    blaze::DynamicVector<double> w = _quad->GetWeights();

    for( int k = 0; k < _nQuadPoints; ++k ) {
        blaze::DynamicVector<double> uKinetic = UKinetic( lambda * _phiTilde[k] );
        for( int j = 0; j < _nMoments; ++j ) {
            for( int l = 0; l < _nStates; ++l ) {
                g[l * _nMoments + j] += uKinetic[l] * _phiTilde[k][j] * w[k];
            }
        }
    }

    return 0.5 * g - MakeVector( u );
}

blaze::DynamicMatrix<double> Closure::Hessian( const blaze::DynamicMatrix<double>& lambda ) {
    blaze::DynamicMatrix<double> H( _nStates * _nMoments, _nStates * _nMoments, 0.0 );
    std::vector<blaze::DynamicMatrix<double>> dLambdadU( _nQuadPoints, blaze::DynamicMatrix<double>( _nStates, _nStates, 0.0 ) );
    for( int k = 0; k < _nQuadPoints; ++k ) {
        dLambdadU[k] = DUKinetic( lambda * _phiTilde[k] );
    }
    for( int l = 0; l < _nStates; ++l ) {
        for( int m = 0; m < _nStates; ++m ) {
            for( int i = 0; i < _nMoments; ++i ) {
                for( int j = 0; j < _nMoments; ++j ) {
                    for( int k = 0; k < _nQuadPoints; ++k ) {    // TODO: reorder to avoid cach misses
                        H( m * _nMoments + j, l * _nMoments + i ) += _hPartial[k]( i, j ) * dLambdadU[k]( l, m );
                    }
                }
            }
        }
    }
    return 0.5 * H;
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
