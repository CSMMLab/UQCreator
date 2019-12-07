#include "thermalpn.h"
#include "quadraturegrid.h"

ThermalPN::ThermalPN( Settings* settings ) : Problem( settings ) {
    _N = 2;
    if( _N == 1 )    // GlobalIndex has different ordering here
        _nMoments = 4;
    else
        _nMoments = unsigned( GlobalIndex( _N, _N ) + 1 );
    _nStates = _nMoments + 1;
    this->SetupSystemMatrices();
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    _settings->SetSource( false );
    _suOlson         = false;
    _constitutiveLaw = 2;    // 1 is Su Olson, 2 is constant

    // physical constants
    _c             = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _a             = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _TRef          = 1.0;                    // reference temperature
    _sigma         = 1.0;                    // opacity
    _alpha         = 4.0 * _a;               // closure relation
    _cV            = 0.718 * 1e7;            // heat capacity. Air has cV = 0.718 in [kJ/(kg K)] = [1000 m^2 / s^2 / K]
    double sigmaSB = 5.6704 * 1e-5;          // Stefan Boltzmann constant in [erg/cm^2/s/K^4]
    _a             = 4.0 * sigmaSB / _c;

    if( !_suOlson ) _TRef = pow( _cV / _a, 1.0 / 4.0 );    // ensure eTilde = O(1)

    _epsilon = 4.0 * _a / _alpha;

    // compute xi Quadrature points
    Vector xiEta( _settings->GetNDimXi() );
    _variances = _settings->GetSigma();

    // get quadrature grid
    auto grid = QuadratureGrid::Create( _settings, _settings->GetNQuadPoints() );
    _xiQuad   = grid->GetNodes();

    // compute Roe flux components
    Matrix P( 2, 2, 1.0 );
    P( 0, 0 ) = -sqrt( 3 ) / _c;
    P( 0, 1 ) = sqrt( 3 ) / _c;
    Matrix PInv( 2, 2, 0.5 );
    PInv( 0, 0 ) = -_c / sqrt( 3 ) / 2.0;
    PInv( 1, 0 ) = _c / sqrt( 3 ) / 2.0;
    Matrix LambdaAbs( 2, 2, 0.0 );
    LambdaAbs( 0, 0 ) = fabs( -1.0 / sqrt( 3 ) / _epsilon );
    LambdaAbs( 1, 1 ) = fabs( 1.0 / sqrt( 3 ) / _epsilon );
    Matrix AbsAPart   = PInv * LambdaAbs * P;
    _AbsA             = Matrix( _nStates, _nStates, 0.0 );
    _AbsA( 0, 0 )     = AbsAPart( 0, 0 );
    _AbsA( 0, 3 )     = AbsAPart( 0, 1 );
    _AbsA( 3, 0 )     = AbsAPart( 1, 0 );
    _AbsA( 3, 3 )     = AbsAPart( 1, 1 );

    // scale refinement threshholds
    //_settings->SetRefinementThreshold( 1e-27 * _settings->GetRefinementThreshold() / ( _a * pow( _TRef, 4 ) ) );
    //_settings->SetCoarsenThreshold( 1e-27 * _settings->GetCoarsenThreshold() / ( _a * pow( _TRef, 4 ) ) );

    _settings->SetRefinementThreshold( 1e-15 * _settings->GetRefinementThreshold() / ( _a * pow( _TRef, 4 ) ) );
    _settings->SetCoarsenThreshold( 1e-15 * _settings->GetCoarsenThreshold() / ( _a * pow( _TRef, 4 ) ) );

    // std::cout << "Threshold " << _settings->GetRefinementThreshold() << " " << _settings->GetCoarsenThreshold() << std::endl;
    // std::cout << _Az << std::endl;
    // exit( EXIT_FAILURE );

    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ThermalPN] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ThermalPN::~ThermalPN() {}

Vector ThermalPN::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    Vector g = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT();
    // Vector g     = 0.5 * ( F( u ) + F( v ) ) * nUnit - 0.5 * _AbsA * ( v - u );
    g[_nMoments] = 0.0;    // set temperature flux to zero
    return g;
}

Matrix ThermalPN::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ThermalPN::F( const Vector& u ) {
    Matrix flux( u.size(), 1, 0.0 );
    Vector outZ = _c / _epsilon * _Az * u;
    outZ[0]     = 1.0 / _c / _c * outZ[0];

    for( unsigned i = 0; i < _nMoments; ++i ) flux( i, 0 ) = outZ[i];
    // flux( _nMoments, 0 ) = 0.0;    // set temperature flux to zero
    return flux;
}

double ThermalPN::Delta( int l, int k ) const {
    if( l == 0 ) {
        return std::sqrt( 4.0 * M_PI );
    }
    else if( l == 1 && k == 0 ) {
        return std::sqrt( 4.0 * M_PI / 3.0 );
    }
    else if( l == 1 && ( k == -1 || k == 1 ) ) {
        return std::sqrt( 4.0 * M_PI * double( MathTools::Factorial( l + std::abs( k ) ) ) /
                          ( 2.0 * ( 2.0 * l + 1 ) * double( MathTools::Factorial( l - std::abs( k ) ) ) ) );
    }
    else {
        return 1.0;
    }
}

void ThermalPN::SetupSystemMatrices() {
    int j;
    unsigned i;
    _Ax = Matrix( _nMoments, _nMoments );
    _Ay = Matrix( _nMoments, _nMoments );
    _Az = Matrix( _nMoments, _nMoments );
    // loop over columns of A
    for( int l = 0; l <= _N; ++l ) {
        for( int k = -l; k <= l; ++k ) {
            i = unsigned( GlobalIndex( l, k ) );

            // flux matrix in direction x
            if( k != -1 ) {
                j = GlobalIndex( l - 1, kMinus( k ) );
                if( j >= 0 && j < int( _nMoments ) )
                    _Ax( i, unsigned( j ) ) = 0.5 * CTilde( l - 1, std::abs( k ) - 1 ) / Delta( l - 1, std::abs( k ) - 1 );
                j = GlobalIndex( l + 1, kMinus( k ) );
                if( j >= 0 && j < int( _nMoments ) ) {
                    _Ax( i, unsigned( j ) ) = -0.5 * DTilde( l + 1, std::abs( k ) - 1 ) / Delta( l + 1, std::abs( k ) - 1 );
                }
            }

            j = GlobalIndex( l - 1, kPlus( k ) );
            if( j >= 0 && j < int( _nMoments ) )
                _Ax( i, unsigned( j ) ) = -0.5 * ETilde( l - 1, std::abs( k ) + 1 ) / Delta( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, kPlus( k ) );
            if( j >= 0 && j < int( _nMoments ) )
                _Ax( i, unsigned( j ) ) = 0.5 * FTilde( l + 1, std::abs( k ) + 1 ) / Delta( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction y
            if( k != 1 ) {
                j = GlobalIndex( l - 1, -kMinus( k ) );
                if( j >= 0 && j < int( _nMoments ) )
                    _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * CTilde( l - 1, std::abs( k ) - 1 ) / Delta( l - 1, std::abs( k ) - 1 );

                j = GlobalIndex( l + 1, -kMinus( k ) );
                if( j >= 0 && j < int( _nMoments ) )
                    _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * DTilde( l + 1, std::abs( k ) - 1 ) / Delta( l + 1, std::abs( k ) - 1 );
            }

            j = GlobalIndex( l - 1, -kPlus( k ) );
            if( j >= 0 && j < int( _nMoments ) )
                _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * ETilde( l - 1, std::abs( k ) + 1 ) / Delta( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, -kPlus( k ) );
            if( j >= 0 && j < int( _nMoments ) )
                _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * FTilde( l + 1, std::abs( k ) + 1 ) / Delta( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction z
            j = GlobalIndex( l - 1, k );
            if( j >= 0 && j < int( _nMoments ) ) _Az( i, unsigned( j ) ) = AParam( l - 1, k ) / Delta( l - 1, k );

            j = GlobalIndex( l + 1, k );
            if( j >= 0 && j < int( _nMoments ) ) _Az( i, unsigned( j ) ) = BParam( l + 1, k ) / Delta( l + 1, k );

            // multiply to change to monomials for up to order one
            for( unsigned n = 0; n < _nMoments; ++n ) {
                _Ax( i, n ) = _Ax( i, n ) * Delta( l, k );
                _Ay( i, n ) = _Ay( i, n ) * Delta( l, k );
                _Az( i, n ) = _Az( i, n ) * Delta( l, k );
            }
        }
    }
}

int ThermalPN::GlobalIndex( int l, int k ) const {
    if( l != 1 ) {
        int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
        int prevIndicesThisLevel = k + l;    // number of previous indices in current level
        return numIndicesPrevLevel + prevIndicesThisLevel;
    }
    else {
        if( k == 1 ) {
            return 1;
        }
        else if( k == -1 ) {
            return 2;
        }
        else {    // k == 0
            return 3;
        }
    }
}

Matrix ThermalPN::Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const {
    unsigned nStates             = static_cast<unsigned>( uQ.rows() );
    unsigned Nq                  = _settings->GetNqPEAtRef( level );
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( level );    // get indices in quadrature array for current refinement level

    Matrix y( nStates, Nq, 0.0 );
    double S           = 0.0;    // source, needs to be defined
    double varianceVal = 0;

    for( unsigned k = 0; k < qIndex.size(); ++k ) {
        if( _suOlson && t < 10 && std::fabs( x[0] ) < 0.5 + _variances[0] * _xiQuad[qIndex[k]][0] ) {
            S           = _a;
            varianceVal = _variances[0];
        }
        else {
            S = 0.0;
        }

        double Q = S / _sigma / _a / std::pow( _TRef, 4 );

        double E      = uQ( 0, k );
        double F      = uQ( 3, k );
        double eTilde = uQ( _nMoments, k );    // scaled internal energy
        if( eTilde < 0 ) {
            std::cout << "eTilde < 0 !!!!" << std::endl;
            std::cout << "eTilde = " << eTilde << std::endl;
            std::cout << "E = " << E << std::endl;
            std::cout << "F = " << F << std::endl;
        }
        double TTilde = ScaledTemperature( eTilde );
        // y( 0, k ) = ( -( E - U ) + ( Q + varianceVal * _xiQuad[k][0] ) ) / _epsilon;

        y( 0, k )         = ( -( E - std::pow( TTilde, 4 ) ) + Q ) / _epsilon;
        y( 3, k )         = -F / _epsilon;
        y( _nMoments, k ) = ( E - std::pow( TTilde, 4 ) ) / _epsilon;
    }

    return y;
}

double ThermalPN::ScaledInternalEnergy( double TTilde ) const {
    double T = TTilde * _TRef;
    // std::cout << "T = " << T << std::endl;
    double e;
    if( _constitutiveLaw == 1 ) {
        e = _alpha / 4.0 * pow( T, 4 );
    }
    else {
        e = _cV * T;
    }
    // std::cout << "cV * T = " << e << std::endl;
    return e / ( _a * pow( _TRef, 4 ) );
}

double ThermalPN::ScaledTemperature( double eTilde ) const {
    double e = eTilde * _a * pow( _TRef, 4 );
    double T;
    if( _constitutiveLaw == 1 ) {
        T = pow( 4.0 * e / _alpha, 1.0 / 4.0 );
    }
    else {
        T = e / _cV;
    }
    return T / _TRef;
}

Matrix ThermalPN::F( const Matrix& u ) {
    _log->error( "[ThermalPN] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ThermalPN::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    double cfl = _settings->GetCFL();

    double maxVelocity = std::sqrt( 1 / 3.0 ) / _epsilon;

    return ( cfl * dx ) / maxVelocity;
}

Vector ThermalPN::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    auto sigma = _settings->GetSigma();
    double E, F, internalEnergy;
    double T = 0;

    if( _suOlson ) {
        E = 0.0;    // std::fmax( 1e-4 * _a,
                    //_a * pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) *
                    //  exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );
        F              = 0;
        internalEnergy = 1e-7 * _a * pow( _TRef, 4 );    // fix to ensure positive values of the inner energy - use 1e-3 without IPM
    }
    else {
        double a = 0.275;
        double b = 0.1;
        if( xi.size() > 1 && false ) {
            a = a + sigma[1] * xi[1];
        }
        if( xi.size() > 2 ) {
            b = b + sigma[2] * xi[2];
        }
        double alpha = pow( a, 1.0 / 4.0 );
        double beta  = pow( b, 1.0 / 4.0 );
        double tau0  = 0.1 + sigma[0] * xi[0];

        if( x[0] < 0.0 )
            T = alpha;
        else if( x[0] < tau0 ) {
            T = 1;
            if( xi.size() > 1 ) {
                T = T + sigma[1] * xi[1];
            }
        }
        else
            T = beta;
        F              = 0.0;
        internalEnergy = ScaledInternalEnergy( T / _TRef ) * ( _a * pow( _TRef, 4 ) );    // 1e-7 * _a * pow( T, 4 );
        E              = std::pow( T / _TRef, 4 );
    }

    y[0] = E / _a / pow( _TRef, 4 );
    if( !_suOlson ) y[0] = std::pow( T / _TRef, 4 );
    y[3]         = F / _a / pow( _TRef, 4 );
    y[_nMoments] = internalEnergy / ( _a * pow( _TRef, 4 ) );
    // std::cout << y[2] << std::endl;
    return y;
}

Vector ThermalPN::LoadIC( const Vector& x, const Vector& xi ) {
    _log->error( "[ThermalPN: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}

double ThermalPN::AParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l + k + 1 ) ) / double( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) );
}

double ThermalPN::BParam( int l, int k ) const { return std::sqrt( double( ( l - k ) * ( l + k ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) ); }

double ThermalPN::CParam( int l, int k ) const {
    return std::sqrt( double( ( l + k + 1 ) * ( l + k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double ThermalPN::DParam( int l, int k ) const {
    return std::sqrt( double( ( l - k ) * ( l - k - 1 ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
}

double ThermalPN::EParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l - k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double ThermalPN::FParam( int l, int k ) const { return std::sqrt( double( ( l + k ) * ( l + k - 1 ) ) / double( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ); }

double ThermalPN::CTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * CParam( l, k );
    else
        return CParam( l, k );
}

double ThermalPN::DTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * DParam( l, k );
    else
        return DParam( l, k );
}

double ThermalPN::ETilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * EParam( l, k );
    else
        return EParam( l, k );
}

double ThermalPN::FTilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * FParam( l, k );
    else {
        return FParam( l, k );
    }
}

int ThermalPN::Sgn( int k ) const {
    if( k >= 0 )
        return 1;
    else
        return -1;
}

int ThermalPN::kPlus( int k ) const { return k + Sgn( k ); }

int ThermalPN::kMinus( int k ) const { return k - Sgn( k ); }
