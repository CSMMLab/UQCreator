#include "closure.h"
#include "closures/boundedbarrier.h"
#include "closures/eulerclosure.h"
#include "closures/eulerclosure2d.h"
#include "closures/hyperbolicitylimiter.h"
#include "closures/hyperbolicitylimiter2d.h"
#include "closures/kineticclosure.h"
#include "closures/logbarrierclosure.h"
#include "closures/logsin.h"
#include "closures/m1ipmclosure.h"
#include "closures/radihydroclosure1d.h"
#include "closures/regularizedeuler1d.h"
#include "closures/regularizedeuler2d.h"
#include "closures/shallowwaterclosure.h"
#include "closures/shallowwaterclosure2d.h"
#include "closures/stochasticgalerkin.h"
#include "closures/thermalradiationclosure.h"
#include "mathtools.h"
#include "tensorizedquadrature.h"
#include "uniformsparsegrid.h"

Closure::Closure( Settings* settings )
    : _settings( settings ), _maxDegree( _settings->GetMaxDegree() ), _nQuadPoints( _settings->GetNQuadPoints() ),
      _nStates( _settings->GetNStates() ), _maxIterations( _settings->GetMaxIterations() ) {
    _log = spdlog::get( "event" );
    // initialize classes: Basis Functions and Quadrature rules are defined for Legendre and Hermite
    _quad.resize( 2 );
    _quad[0] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_LEGENDRE );
    _quad[1] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_HERMITE );

    // get number of uncertain dimensions
    _numDimXi = _settings->GetNDimXi();

    // get number of multi elements per cell
    _nMultiElements = _settings->GetNMultiElements();

    // get index vector for basis function calculation
    std::vector<std::vector<unsigned>> indices = _settings->GetPolyIndices();
    _nTotal                                    = _settings->GetNTotal();

    // define quadrature
    _wGrid.resize( _settings->GetNRefinementLevels() );
    _nQTotalForRef.resize( _settings->GetNRefinementLevels() );
    _nTotalForRef      = _settings->GetNTotalRefinementLevel();
    VectorU quadLevel  = _settings->GetQuadLevel();
    unsigned oldQLevel = 0;
    _quadGrid          = nullptr;
    for( unsigned rlevel = 0; rlevel < _settings->GetNRefinementLevels(); ++rlevel ) {
        if( oldQLevel != quadLevel[rlevel] ) {    // new quadrature weights needed
            if( _quadGrid ) delete _quadGrid;
            _quadGrid = QuadratureGrid::Create( _settings, quadLevel[rlevel] );
        }
        _wGrid[rlevel]         = _quadGrid->GetWeights();
        _nQTotalForRef[rlevel] = _wGrid[rlevel].size();
        oldQLevel              = quadLevel[rlevel];
    }
    _xiGrid = _quadGrid->GetNodes();
    // give number of quad points per level to settings. Quadrature variables are computed here.
    _settings->SetNQTotalForRef( _nQTotalForRef );

    // set total number of quadrature points
    _nQTotal = _quadGrid->GetNodeCount();

    // compute basis functions evaluated at the quadrature points
    Matrix phiTilde = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeF      = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeWf     = Matrix( _nQTotal, _nTotal, 1.0 );
    _phiTildeVec    = std::vector<Vector>( _nQTotal, Vector( _nTotal, 0.0 ) );

    unsigned n = 0;
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned m = 0; m < _numDimXi; ++m ) {
                if( _settings->GetDistributionType( m ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( m ) == DistributionType::D_HERMITE ) n = 1;
                phiTilde( k, i ) *=
                    _quad[n]->Evaluate( indices[i][m], _xiGrid[k][m] ) / _quad[n]->L2Norm( indices[i][m] );    // sqrt( 2.0 * i + 1.0 );
                _phiTildeF( k, i ) *=
                    _quad[n]->Evaluate( indices[i][m], _xiGrid[k][m] ) / _quad[n]->L2Norm( indices[i][m] ) * _quad[n]->fXi( _xiGrid[k][m] );
            }
            // P(xi \in I_l) probability that xi lies in multi element I_l
            // double P            = 1.0 / _nMultiElements;
            //_phiTildeF( k, i )  = _phiTildeF( k, i ) / P;    // modify pdf to multielement ansatz, only valid for uniform distributions
            _phiTildeWf( k, i ) = _phiTildeF( k, i ) * _wGrid[_settings->GetNRefinementLevels() - 1][k];
            _phiTildeVec[k][i]  = phiTilde( k, i );
        }
    }

    _phiTildeTrans      = trans( phiTilde );
    auto phiTildeFTrans = trans( _phiTildeF );
    // calculate partial matrix for Hessian calculation
    _hPartial = MatVec( _nQTotal, Matrix( _nTotal, _nTotal ) );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        _hPartial[k] = outer( column( _phiTildeTrans, k ), column( phiTildeFTrans, k ) );    // TODO
    }

    _filter = Filter::Create( _settings );
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
        std::cout << "test I " << testQuad << std::endl;*/
}

Closure::~Closure() {
    for( unsigned l = 0; l < _quad.size(); ++l ) {
        delete _quad[l];
    }
    delete _quadGrid;
    delete _filter;
}

Closure* Closure::Create( Settings* settings ) {
    auto log         = spdlog::get( "event" );
    auto closureType = settings->GetClosureType();
    if( closureType == ClosureType::C_BOUNDEDBARRIER ) {
        return new BoundedBarrier( settings );
    }
    else if( closureType == ClosureType::C_LOGBARRIER ) {
        return new LogBarrierClosure( settings );
    }
    else if( closureType == ClosureType::C_LOGSIN ) {
        return new LogSin( settings );
    }
    else if( closureType == ClosureType::C_STOCHASTICGALERKIN ) {
        return new StochasticGalerkin( settings );
    }
    else if( closureType == ClosureType::C_EULER_1D ) {
        if( settings->HasRegularization() ) {
            return new RegularizedEuler1D( settings );
        }
        else {
            return new EulerClosure( settings );
        }
    }
    else if( closureType == ClosureType::C_EULER_2D ) {
        if( settings->HasRegularization() ) {
            return new RegularizedEuler2D( settings );
        }
        else {
            return new EulerClosure2D( settings );    // no filtering, no regularization
        }
    }
    else if( closureType == ClosureType::C_SHALLOWWATER_1D ) {
        return new ShallowWaterClosure( settings );
    }
    else if( closureType == ClosureType::C_SHALLOWWATER_2D ) {
        return new ShallowWaterClosure2D( settings );
    }
    else if( closureType == ClosureType::C_RADHYDRO ) {
        return new RadiHydroClosure1D( settings );
    }
    else if( closureType == ClosureType::C_THERMALRAD_1D ) {
        return new ThermalRadiationClosure( settings );
    }
    else if( closureType == ClosureType::C_M1_1D ) {
        return new M1IPMClosure( settings );
    }
    else if( closureType == ClosureType::C_KINETIC ) {
        return new KineticClosure( settings );
    }
    else if( closureType == ClosureType::C_HYPLIM ) {
        return new HyperbolicityLimiter( settings );
    }
    else if( closureType == ClosureType::C_HYPLIM_2D ) {
        return new HyperbolicityLimiter2D( settings );
    }
    else {
        log->error( "[closure]: Invalid closure type" );
        exit( EXIT_FAILURE );
    }
}

void Closure::SolveClosure( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
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
        Gradient( g, lambdaMat, uMat, refLevel );
        if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
            continue;
        }

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

void Closure::SolveClosureSafe( Tensor& lambda, const Tensor& u, unsigned refLevel ) {
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
        Gradient( g, lambdaMat, uMat, refLevel );
        if( CalcNorm( g, nTotal ) < _settings->GetEpsilon() ) {
            continue;
        }

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
        // std::cout << CalcNorm( dlambdaNew, nTotal ) << std::endl;
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
                // std::cout << "Refining: " << CalcNorm( dlambdaNew, nTotal ) << std::endl;
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

double Closure::CalcNorm( Vector& test, unsigned nTotal ) const {
    double out = 0.0;
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out += pow( test[s * nTotal + i], 2 );
        }
    }
    return sqrt( out );
}

Vector Closure::EvaluateLambda( const Tensor& lambda, unsigned l, unsigned k, unsigned nTotal ) {
    Vector out( _nStates, 0.0 );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out[s] += lambda( s, l, i ) * _phiTildeVec[k][i];
        }
    }
    return out;
}

Vector Closure::EvaluateLambda( const Matrix& lambda, unsigned k, unsigned nTotal ) {
    Vector out( _nStates, 0.0 );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out[s] += lambda( s, i ) * _phiTildeVec[k][i];
        }
    }
    return out;
}

Tensor Closure::EvaluateLambda( const Tensor& lambda ) const { return lambda * _phiTildeTrans; }

Tensor Closure::EvaluateLambdaOnPE( const Tensor& lambda, unsigned levelOld, unsigned levelNew ) const {
    Tensor out( _settings->GetNStates(), _nMultiElements, _settings->GetNqPEAtRef( levelNew ), 0.0 );
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( levelNew );
    unsigned nTotal              = _nTotalForRef[levelOld];

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned l = 0; l < _nMultiElements; ++l ) {
            for( unsigned k = 0; k < qIndex.size(); ++k ) {
                for( unsigned i = 0; i < nTotal; ++i ) {
                    out( s, l, k ) += lambda( s, l, i ) * _phiTildeTrans( i, qIndex[k] );
                }
            }
        }
    }
    return out;
}

void Closure::EvaluateLambda( Tensor& out, const Tensor& lambda ) const { out = lambda * _phiTildeTrans; }

void Closure::Gradient( Vector& g, const Matrix& lambda, const Matrix& u, unsigned refLevel ) {
    Vector uKinetic( _nStates, 0.0 );
    unsigned nTotal = _nTotalForRef[refLevel];
    g.reset();

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {
        U( uKinetic, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned i = 0; i < nTotal; ++i ) {
            for( unsigned s = 0; s < _nStates; ++s ) {
                g[s * nTotal + i] += uKinetic[s] * _phiTildeF( k, i ) * _wGrid[refLevel][k];
            }
        }
    }

    SubstractVectorMatrixOnVector( g, u, _nTotalForRef[refLevel] );
}

void Closure::Hessian( Matrix& H, const Matrix& lambda, unsigned refLevel ) {
    H.reset();
    Matrix dUdLambda( _nStates, _nStates );    // TODO: preallocate Matrix for Hessian computation -> problems omp
    unsigned nTotal = _nTotalForRef[refLevel];

    for( unsigned k = 0; k < _nQTotalForRef[refLevel]; ++k ) {    // TODO: reorder to avoid cache misses
        DU( dUdLambda, EvaluateLambda( lambda, k, nTotal ) );
        for( unsigned s = 0; s < _nStates; ++s ) {
            for( unsigned m = 0; m < _nStates; ++m ) {
                for( unsigned j = 0; j < nTotal; ++j ) {
                    for( unsigned i = 0; i < nTotal; ++i ) {
                        H( m * nTotal + j, s * nTotal + i ) += _hPartial[k]( j, i ) * _wGrid[refLevel][k] * dUdLambda( s, m );
                    }
                }
            }
        }
    }
}

void Closure::AddMatrixVectorToMatrix( const Matrix& A, const Vector& b, Matrix& y, unsigned nTotal ) const {
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            y( s, j ) = A( s, j ) + b[s * nTotal + j];
        }
    }
}

void Closure::SubstractVectorMatrixOnVector( Vector& b, const Matrix& A, unsigned nTotal ) const {
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned j = 0; j < nTotal; ++j ) {
            b[s * nTotal + j] = b[s * nTotal + j] - A( s, j );
        }
    }
}

void Closure::DS( Vector& /*ds*/, const Vector& /*u*/ ) const {}

std::vector<Polynomial*> Closure::GetQuadrature() { return _quad; }

void Closure::SetAlpha( double alpha ) { _alpha = alpha; }

void Closure::SetMaxIterations( unsigned maxIterations ) { _maxIterations = maxIterations; }

unsigned Closure::GetMaxIterations() const { return _maxIterations; }

QuadratureGrid* Closure::GetQuadratureGrid() { return _quadGrid; }

Matrix Closure::GetPhiTildeWfAtRef( unsigned level ) const {
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( level );
    Matrix phiTildeWfTrans( unsigned( qIndex.size() ), _nTotalForRef[level], false );
    for( unsigned k = 0; k < qIndex.size(); ++k ) {
        for( unsigned i = 0; i < _nTotalForRef[level]; ++i ) {
            phiTildeWfTrans( k, i ) = _phiTildeF( qIndex[k], i ) * _wGrid[level][qIndex[k]];
        }
    }
    return phiTildeWfTrans;
}

Matrix Closure::GetPhiTildeWfAtRef( unsigned level, bool /*full*/ ) const {
    unsigned kStart = 0;
    unsigned kEnd   = _nQTotalForRef[level] - 1;
    Matrix phiTildeWfTrans( _nQTotalForRef[level], _nTotalForRef[level], false );
    for( unsigned k = kStart; k <= kEnd; ++k ) {
        for( unsigned i = 0; i < _nTotalForRef[level]; ++i ) {
            phiTildeWfTrans( k - kStart, i ) = _phiTildeF( k, i ) * _wGrid[level][k];
        }
    }
    return phiTildeWfTrans;
}
