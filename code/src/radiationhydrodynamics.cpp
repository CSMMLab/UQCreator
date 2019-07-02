#include "radiationhydrodynamics.h"

RadiationHydrodynamics::RadiationHydrodynamics( Settings* settings ) : PNEquations( settings, true ) {
    _N        = 1;    // set moment order
    _c        = 1.0;
    _P        = 1.0;
    _gamma    = 1.4;
    _R        = 287.87;
    _nMoments = unsigned( GlobalIndex( _N, _N ) + 1 );
    _nStates  = _nMoments + 4;    // total number of states in equations in P_N + Euler
    _settings->SetNStates( _nStates );
    _settings->SetSource( false );    // TODO:DEBUGGING
    _sigmaA = 0.0;                    // absorption coefficient
    _sigmaS = 0.0;                    // scattering coefficient
    _sigmaT = _sigmaA + _sigmaS;
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
        SetupSystemMatrices();
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[pnequations] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Matrix RadiationHydrodynamics::Source( const Matrix& uQ ) const {
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = static_cast<unsigned>( uQ.columns() );
    Matrix y( nStates, Nq, 0.0 );
    for( unsigned k = 0; k < Nq; ++k ) {
        double se             = SE( column( uQ, k ) );
        Vector sFVec          = SF( column( uQ, k ) );
        y( 0, k )             = _c * se;
        y( 1, k )             = _c * sFVec[0];
        y( 2, k )             = _c * sFVec[1];
        y( 3, k )             = _c * sFVec[2];
        y( _nMoments + 1, k ) = -_P * sFVec[0];
        y( _nMoments + 2, k ) = -_P * sFVec[1];
        y( _nMoments + 3, k ) = -_P * _c * se;
    }
    return y;
}

double RadiationHydrodynamics::Delta( int l, int k ) const {
    return 1.0;    // TODO: DEBUGGING
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

void RadiationHydrodynamics::SetupSystemMatrices() {
    int j;
    unsigned i;
    unsigned nTotalEntries = unsigned( GlobalIndex( _N, _N ) + 1 );    // total number of entries for sytem matrix
    _Ax                    = Matrix( nTotalEntries, nTotalEntries );
    _Ay                    = Matrix( nTotalEntries, nTotalEntries );
    _Az                    = Matrix( nTotalEntries, nTotalEntries );
    // loop over columns of A
    for( int l = 0; l <= _N; ++l ) {
        for( int k = -l; k <= l; ++k ) {
            i = unsigned( GlobalIndex( l, k ) );

            // flux matrix in direction x
            if( k != -1 ) {
                j = GlobalIndex( l - 1, kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) )
                    _Ax( i, unsigned( j ) ) = 0.5 * Delta( l - 1, std::abs( k ) - 1 ) * CTilde( l - 1, std::abs( k ) - 1 );
                j = GlobalIndex( l + 1, kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) ) {
                    _Ax( i, unsigned( j ) ) = -0.5 * Delta( l + 1, std::abs( k ) - 1 ) * DTilde( l + 1, std::abs( k ) - 1 );
                }
            }

            j = GlobalIndex( l - 1, kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) )
                _Ax( i, unsigned( j ) ) = -0.5 * Delta( l - 1, std::abs( k ) + 1 ) * ETilde( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) )
                _Ax( i, unsigned( j ) ) = 0.5 * Delta( l + 1, std::abs( k ) + 1 ) * FTilde( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction y
            if( k != 1 ) {
                j = GlobalIndex( l - 1, -kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) )
                    _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * Delta( l - 1, std::abs( k ) - 1 ) * CTilde( l - 1, std::abs( k ) - 1 );

                j = GlobalIndex( l + 1, -kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) )
                    _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * Delta( l + 1, std::abs( k ) - 1 ) * DTilde( l + 1, std::abs( k ) - 1 );
            }

            j = GlobalIndex( l - 1, -kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) )
                _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * Delta( l - 1, std::abs( k ) + 1 ) * ETilde( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, -kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) )
                _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * Delta( l + 1, std::abs( k ) + 1 ) * FTilde( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction z
            j = GlobalIndex( l - 1, k );
            if( j >= 0 && j < int( nTotalEntries ) ) _Az( i, unsigned( j ) ) = Delta( l - 1, k ) * AParam( l - 1, k );

            j = GlobalIndex( l + 1, k );
            if( j >= 0 && j < int( nTotalEntries ) ) _Az( i, unsigned( j ) ) = Delta( l + 1, k ) * BParam( l + 1, k );

            // multiply to change to monomials for up to order one
            for( unsigned n = 0; n < nTotalEntries; ++n ) {
                _Ax( i, n ) = _c * _Ax( i, n ) * Delta( l, k );
                _Ay( i, n ) = _c * _Ay( i, n ) * Delta( l, k );
                _Az( i, n ) = _c * _Az( i, n ) * Delta( l, k );
            }
        }
    }
}
/*
int RadiationHydrodynamics::GlobalIndex( int l, int k ) const {
    if( l != 1 ) {
        int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
        int prevIndicesThisLevel = k + l;    // number of previous indices in current level
        return numIndicesPrevLevel + prevIndicesThisLevel;
    }
    else {
        if( k == 0 ) {
            return 1;
        }
        else if( k == -1 ) {
            return 2;
        }
    }
}*/

double RadiationHydrodynamics::Er0( const Vector& u ) const {
    double Er     = u[0];
    double rhoInv = 1.0 / u[_nMoments + 0];
    double uU     = u[_nMoments + 1] * rhoInv;
    double vU     = u[_nMoments + 2] * rhoInv;
    double wU     = 0.0;
    return Er - 2.0 * ( uU * u[1] + vU * u[2] + wU * u[3] ) / _c;
}

Vector RadiationHydrodynamics::Fr0( const Vector& u ) const {
    double rhoInv = 1.0 / u[_nMoments + 0];
    double Er     = u[0];
    Vector Fr( 3, false );
    Vector v( 3, false );
    Fr[0] = u[1];
    Fr[1] = u[2];
    Fr[2] = u[3];
    v[0]  = u[_nMoments + 1] * rhoInv;
    v[1]  = u[_nMoments + 2] * rhoInv;
    v[2]  = 0.0;
    return Fr - ( v * Er + _Ax * u + _Ay * u + _Az * u ) / _c;
}

double RadiationHydrodynamics::SE( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0]     = u[_nMoments + 1] / rho;
    v[1]     = u[_nMoments + 2] / rho;
    v[2]     = 0.0;
    double p = ( _gamma - 1.0 ) * ( u[_nMoments + 3] - 0.5 * rho * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    double T = p / ( _R * rho );
    return _sigmaA * ( std::pow( T, 4 ) - Er ) + ( _sigmaA - _sigmaS ) * v.inner( Fr0( u ) ) / _c;
}

Vector RadiationHydrodynamics::SF( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0]     = u[_nMoments + 1] / rho;
    v[1]     = u[_nMoments + 2] / rho;
    v[2]     = 0.0;
    double p = ( _gamma - 1.0 ) * ( u[_nMoments + 3] - 0.5 * rho * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    double T = p / ( _R * rho );
    return Fr0( u ) * ( -_sigmaT ) + v * _sigmaA * ( std::pow( T, 4 ) - Er ) / _c;
}

Matrix RadiationHydrodynamics::F( const Vector& u ) {
    Matrix flux( u.size(), 2 );
    double rhoInv = 1.0 / u[_nMoments + 0];
    double v1     = u[_nMoments + 1] * rhoInv;
    double v2     = u[_nMoments + 2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[_nMoments + 3] - 0.5 * u[_nMoments + 0] * ( pow( v1, 2 ) + pow( v2, 2 ) ) );

    Vector momentFluxX = _Ax * u;
    Vector momentFluxY = _Ay * u;
    for( unsigned i = 0; i < _nMoments; ++i ) {
        flux( i, 0 ) = momentFluxX[i];
        flux( i, 1 ) = momentFluxY[i];
    }
    flux( _nMoments + 0, 0 ) = u[_nMoments + 1];
    flux( _nMoments + 1, 0 ) = u[_nMoments + 1] * v1 + p;
    flux( _nMoments + 2, 0 ) = u[_nMoments + 1] * v2;
    flux( _nMoments + 3, 0 ) = ( u[_nMoments + 3] + p ) * v1;
    flux( _nMoments + 0, 1 ) = u[_nMoments + 2];
    flux( _nMoments + 1, 1 ) = u[_nMoments + 2] * v1;
    flux( _nMoments + 2, 1 ) = u[_nMoments + 2] * v2 + p;
    flux( _nMoments + 3, 1 ) = ( u[_nMoments + 3] + p ) * v2;
    return flux;
}

Vector RadiationHydrodynamics::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double x0    = 0.0;
    double y0    = 0.0;
    double s2    = 3.2 * std::pow( 0.01, 2 );    // std::pow( 0.03, 2 );
    double floor = 0.0;
    _sigma       = _settings->GetSigma();

    y[0] = Delta( 0, 0 ) *
           std::fmax( floor, 1.0 / ( 4.0 * M_PI * s2 ) * exp( -( ( x[0] - x0 ) * ( x[0] - x0 ) + ( x[1] - y0 ) * ( x[1] - y0 ) ) / 4.0 / s2 ) );
    y[_nMoments + 0] = 1.0;
    y[_nMoments + 3] = 1.0;
    return y;
}
/*
void RadiationHydrodynamics::DS( Vector& ds, const Vector& u ) const {
    double gamma      = _gamma;
    double rho        = u[_nMoments + 0];
    double rhoU       = u[_nMoments + 1];
    double rhoV       = u[_nMoments + 2];
    double rhoV2      = pow( rhoV, 2 );
    double rhoU2      = pow( rhoU, 2 );
    double rhoE       = u[_nMoments + 3];
    ds[_nMoments + 0] = ( rhoU2 + rhoV2 + gamma * ( 2 * rho * rhoE - rhoU2 - rhoV2 ) ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) -
                        std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );
    ds[_nMoments + 1] = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[_nMoments + 2] = -( ( 2 * rho * rhoV ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[_nMoments + 3] = -( rho / ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );

    // safety factor on ds3
    // ds[3] -= 1e-7;
}*/

Matrix RadiationHydrodynamics::FRadiation( const Vector& u ) {
    Matrix flux( u.size(), 2 );

    column( flux, 0 ) = _Ax * u;
    column( flux, 1 ) = _Ay * u;
    // column( flux, 1 ) = _Az * u;

    return flux;
}

Matrix RadiationHydrodynamics::FEuler( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v1     = u[1] * rhoInv;
    double v2     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( v1, 2 ) + pow( v2, 2 ) ) );
    Matrix flux( u.size(), 2 );

    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v1 + p;
    flux( 2, 0 ) = u[1] * v2;
    flux( 3, 0 ) = ( u[3] + p ) * v1;
    flux( 0, 1 ) = u[2];
    flux( 1, 1 ) = u[2] * v1;
    flux( 2, 1 ) = u[2] * v2 + p;
    flux( 3, 1 ) = ( u[3] + p ) * v2;

    return flux;
}

Vector RadiationHydrodynamics::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    Vector out( _nStates, false );
    Vector outEuler( 4, false );
    Vector uRadiation( _nStates, false );
    Vector vRadiation( _nStates, false );
    // copy radiation part on uRadiation
    for( unsigned i = 0; i < _nMoments; ++i ) {
        uRadiation[i] = u[i];
        vRadiation[i] = v[i];
    }
    double rhoInv = 1.0 / u[_nMoments + 0];
    double uU     = u[_nMoments + 1] * rhoInv;
    double vU     = u[_nMoments + 2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[_nMoments + 3] - 0.5 * u[_nMoments + 0] * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
    double aU     = sqrt( _gamma * p * rhoInv );

    rhoInv    = 1.0 / v[_nMoments + 0];
    double uV = v[_nMoments + 1] * rhoInv;
    double vV = v[_nMoments + 2] * rhoInv;
    p         = ( _gamma - 1.0 ) * ( v[_nMoments + 3] - 0.5 * v[_nMoments + 0] * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
    double aV = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
    double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        outEuler = FEuler( u ) * n;
    else if( lambdaMax <= 0 )
        outEuler = FEuler( v ) * n;
    else {
        outEuler = ( 1.0 / ( lambdaMax - lambdaMin ) ) *
                   ( lambdaMax * FEuler( u ) * n - lambdaMin * FEuler( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }

    // write radiation part on _nMoments entries
    out = FRadiation( 0.5 * ( uRadiation + vRadiation ) ) * n - 0.5 * ( vRadiation - uRadiation ) * norm( n );

    // save Euler part on return vector
    for( unsigned s = 0; s < 4; ++s ) out[_nMoments + s] = outEuler[s];    // TODO: DEBUG

    return out;
}

Matrix RadiationHydrodynamics::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}
