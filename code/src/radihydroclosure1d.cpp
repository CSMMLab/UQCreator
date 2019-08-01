#include "radihydroclosure1d.h"

RadiHydroClosure1D::RadiHydroClosure1D( Settings* settings ) : Closure( settings ), _gamma( _settings->GetGamma() ) {
    _alpha        = 0.5;
    _nMoments     = 4;
    _nHydroStates = 3;
    _gamma        = 5.0 / 3.0;
}

RadiHydroClosure1D::~RadiHydroClosure1D() {}

void RadiHydroClosure1D::U( Vector& out, const Vector& Lambda ) {
    for( unsigned s = 0; s < _nMoments; ++s ) out[s] = Lambda[s];
    out[_nMoments + 0] = exp( ( 2.0 * Lambda[_nMoments + 0] * Lambda[_nMoments + 2] - 2.0 * Lambda[_nMoments + 2] * log( -Lambda[_nMoments + 2] ) -
                                2.0 * Lambda[_nMoments + 2] * _gamma - pow( Lambda[_nMoments + 1], 2.0 ) ) /
                              ( 2.0 * Lambda[_nMoments + 2] * ( _gamma - 1.0 ) ) );
    out[_nMoments + 1] = -( Lambda[_nMoments + 1] / Lambda[_nMoments + 2] ) * out[_nMoments + 0];
    out[_nMoments + 2] =
        ( ( pow( Lambda[_nMoments + 1], 2.0 ) - 2.0 * Lambda[_nMoments + 2] ) / ( 2.0 * pow( Lambda[_nMoments + 2], 2.0 ) ) ) * out[_nMoments + 0];
}

void RadiHydroClosure1D::U( Vector& out, const Vector& Lambda, bool dummy ) {
    out[0] = exp( ( 2.0 * Lambda[0] * Lambda[2] - 2.0 * Lambda[2] * log( -Lambda[2] ) - 2.0 * Lambda[2] * _gamma - pow( Lambda[1], 2.0 ) ) /
                  ( 2.0 * Lambda[2] * ( _gamma - 1.0 ) ) );
    out[1] = -( Lambda[1] / Lambda[2] ) * out[0];
    out[2] = ( ( pow( Lambda[1], 2.0 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2.0 ) ) ) * out[0];
}

void RadiHydroClosure1D::U( Matrix& out, const Matrix& Lambda ) {
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        for( unsigned s = 0; s < _nMoments; ++s ) out( s, k ) = Lambda( s, k );
        out( _nMoments + 0, k ) = exp( ( 2.0 * Lambda( _nMoments + 0, k ) * Lambda( _nMoments + 2, k ) -
                                         2.0 * Lambda( _nMoments + 2, k ) * log( -Lambda( _nMoments + 2, k ) ) -
                                         2.0 * Lambda( _nMoments + 2, k ) * _gamma - pow( Lambda( _nMoments + 1, k ), 2.0 ) ) /
                                       ( 2.0 * Lambda( _nMoments + 2, k ) * ( _gamma - 1.0 ) ) );
        out( _nMoments + 1, k ) = -( Lambda( _nMoments + 1, k ) / Lambda( _nMoments + 2, k ) ) * out( _nMoments + 0, k );
        out( _nMoments + 2, k ) =
            ( ( pow( Lambda( _nMoments + 1, k ), 2.0 ) - 2.0 * Lambda( _nMoments + 2, k ) ) / ( 2.0 * pow( Lambda( _nMoments + 2, k ), 2.0 ) ) ) *
            out( _nMoments + 0, k );
    }
}

Matrix RadiHydroClosure1D::U( const Matrix& Lambda ) {
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        for( unsigned s = 0; s < _nMoments; ++s ) y( s, k ) = Lambda( s, k );
        y( _nMoments + 0, k ) = exp( ( 2.0 * Lambda( _nMoments + 0, k ) * Lambda( _nMoments + 2, k ) -
                                       2.0 * Lambda( _nMoments + 2, k ) * log( -Lambda( _nMoments + 2, k ) ) -
                                       2.0 * Lambda( _nMoments + 2, k ) * _gamma - pow( Lambda( _nMoments + 1, k ), 2.0 ) ) /
                                     ( 2.0 * Lambda( _nMoments + 2, k ) * ( _gamma - 1.0 ) ) );
        y( _nMoments + 1, k ) = -( Lambda( _nMoments + 1, k ) / Lambda( _nMoments + 2, k ) ) * y( _nMoments + 0, k );
        y( _nMoments + 2, k ) =
            ( ( pow( Lambda( _nMoments + 1, k ), 2.0 ) - 2.0 * Lambda( _nMoments + 2, k ) ) / ( 2.0 * pow( Lambda( _nMoments + 2, k ), 2.0 ) ) ) *
            y( _nMoments + 0, k );
    }

    return y;
}

void RadiHydroClosure1D::DU( Matrix& y, const Vector& Lambda ) {
    double E     = exp( ( 2.0 * Lambda[0] * Lambda[2] - 2.0 * Lambda[2] * log( -Lambda[2] ) - 2.0 * Lambda[2] * _gamma - pow( Lambda[1], 2 ) ) /
                    ( 2.0 * Lambda[2] * ( _gamma - 1.0 ) ) );
    double dEdv1 = ( 1.0 / ( _gamma - 1.0 ) ) * E;
    double dEdv2 = -( Lambda[1] / ( Lambda[2] * ( _gamma - 1.0 ) ) ) * E;
    double dEdv3 = ( ( -2.0 * _gamma + 2.0 * Lambda[0] - 2.0 * log( -Lambda[2] ) - 2.0 ) / ( 2.0 * ( _gamma - 1.0 ) * Lambda[2] ) -
                     ( -2.0 * _gamma * Lambda[2] + 2.0 * Lambda[0] * Lambda[2] - pow( Lambda[1], 2 ) - 2.0 * Lambda[2] * log( -Lambda[2] ) ) /
                         ( 2.0 * ( _gamma - 1.0 ) * pow( Lambda[2], 2 ) ) ) *
                   E;
    double preFac3 = ( ( pow( Lambda[1], 2 ) - 2.0 * Lambda[2] ) / ( 2.0 * pow( Lambda[2], 2 ) ) );
    y( 0, 0 )      = dEdv1;
    y( 0, 1 )      = dEdv2;
    y( 0, 2 )      = dEdv3;
    y( 1, 0 )      = -( Lambda[1] / Lambda[2] ) * dEdv1;
    y( 1, 1 )      = -( 1.0 / Lambda[2] ) * E - ( Lambda[1] / Lambda[2] ) * dEdv2;
    y( 1, 2 )      = ( Lambda[1] / ( pow( Lambda[2], 2 ) ) ) * E - ( Lambda[1] / Lambda[2] ) * dEdv3;
    y( 2, 0 )      = preFac3 * dEdv1;
    y( 2, 1 )      = preFac3 * dEdv2 + ( 2.0 * Lambda[1] / ( 2.0 * pow( Lambda[2], 2 ) ) ) * E;
    y( 2, 2 )      = ( -pow( Lambda[1], 2 ) * pow( Lambda[2], -3 ) + pow( Lambda[2], -2 ) ) * E + preFac3 * dEdv3;
}

void RadiHydroClosure1D::DS( Vector& ds, const Vector& u ) const {
    double gamma      = -_gamma;
    double rho        = u[_nMoments + 0];
    double rhoU       = u[_nMoments + 1];
    double rhoU2      = pow( rhoU, 2 );
    double rhoE       = u[_nMoments + 2];
    ds[_nMoments + 0] = ( rhoU2 + gamma * ( 2 * rho * rhoE - rhoU2 ) ) / ( -2 * rho * rhoE + rhoU2 ) -
                        std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 ) / ( 2 * rho ) ) );
    ds[_nMoments + 1] = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 ) );
    ds[_nMoments + 2] = -( rho / ( rhoE - ( rhoU2 ) / ( 2 * rho ) ) );
    for( unsigned s = 0; s < _nMoments; ++s ) ds[s] = u[s];
}

void RadiHydroClosure1D::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
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

void RadiHydroClosure1D::Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel ) {
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

void RadiHydroClosure1D::SolveClosure( Matrix& lambdaFull, const Matrix& uFull, unsigned refLevel ) {
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
        while( CalcNorm( dlambda, nTotal ) < CalcNorm( dlambdaNew, nTotal ) ) {
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

void RadiHydroClosure1D::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            y( l, j ) = A( l, j ) + b[l * nTotal + j];
        }
    }
}

void RadiHydroClosure1D::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const {
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            b[l * nTotal + j] = b[l * nTotal + j] - A( l, j );
        }
    }
}

double RadiHydroClosure1D::CalcNorm( Vector& test, unsigned nTotal ) const {
    double out = 0.0;
    for( unsigned l = 0; l < _nHydroStates; ++l ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out += pow( test[l * nTotal + i], 2 );
        }
    }
    return sqrt( out );
}
