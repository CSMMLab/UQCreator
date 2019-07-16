#include "radiationhydrodynamics1d.h"

RadiationHydrodynamics1D::RadiationHydrodynamics1D( Settings* settings ) : PNEquations( settings, true ) {
    _N            = 1;    // set moment order
    _gamma        = 5.0 / 3.0;
    _R            = 287.87;
    double cLight = 299792458.0 * 100.0;    // speed of light in [cm/s]
    // initial chock states
    double lRef = 1000;
    double rhoL = 5.45887 * 1e-13;
    double uL   = 2.3545 * 1e5;
    double TL   = 100.0;
    double rhoR = 1.2479 * 1e-12;
    double uR   = 1.03 * 1e5;
    double TR   = 207.757;
    double pRef = TL * ( _R * rhoL );
    double aRef = sqrt( _gamma * pRef / rhoL );
    _c          = cLight / aRef;
    _P          = 0.0;    // 0.5;

    if( _N == 1 )    // GlobalIndex has different ordering here
        _nMoments = 4;
    else
        _nMoments = unsigned( GlobalIndex( _N, _N ) + 1 );
    _nStates = _nMoments + 3;    // total number of states in equations in P_N + Euler
    _settings->SetNStates( _nStates );
    _settings->SetSource( true );    // TODO:DEBUGGING
    _sigmaA = 3.93 * 1e-5 * lRef;    // absorption coefficient
    _sigmaT = 0.848902 * lRef;       // 7.5657 * 1e-15;
    _sigmaS = _sigmaT - _sigmaA;     // scattering coefficient
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
        SetupSystemMatrices();
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[pnequations] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Matrix RadiationHydrodynamics1D::Source( const Matrix& uQ ) const {
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
        y( _nMoments + 2, k ) = -_P * _c * se;
    }
    return y;
}

double RadiationHydrodynamics1D::Delta( int l, int k ) const {
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

void RadiationHydrodynamics1D::SetupSystemMatrices() {
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
                _Ax( i, n ) = _c * _Ax( i, n ) * Delta( l, k );
                _Ay( i, n ) = _c * _Ay( i, n ) * Delta( l, k );
                _Az( i, n ) = _c * _Az( i, n ) * Delta( l, k );
            }
        }
    }
}

int RadiationHydrodynamics1D::GlobalIndex( int l, int k ) const {
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

double RadiationHydrodynamics1D::Er0( const Vector& u ) const {
    double Er     = u[0];
    double rhoInv = 1.0 / u[_nMoments + 0];
    double uU     = u[_nMoments + 1] * rhoInv;
    double vU     = 0.0;
    double wU     = 0.0;
    return Er - 2.0 * ( uU * u[1] + vU * u[2] + wU * u[3] ) / _c;
}

Vector RadiationHydrodynamics1D::Fr0( const Vector& u ) const {
    double rhoInv = 1.0 / u[_nMoments + 0];
    double Er     = u[0];
    Vector Fr( 3, false );
    Vector v( 3, false );
    Fr[0] = u[1];
    Fr[1] = u[2];
    Fr[2] = u[3];
    v[0]  = u[_nMoments + 1] * rhoInv;
    v[1]  = 0.0;
    v[2]  = 0.0;
    return Fr - ( v * Er + _Ax * u + _Ay * u + _Az * u ) / _c;
}

double RadiationHydrodynamics1D::SE( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0]     = u[_nMoments + 1] / rho;
    v[1]     = 0.0;
    v[2]     = 0.0;
    double p = ( _gamma - 1.0 ) * ( u[_nMoments + 2] - 0.5 * rho * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    double T = p / ( _R * rho );
    return _sigmaA * ( std::pow( T, 4 ) - Er ) + ( _sigmaA - _sigmaS ) * v.inner( Fr0( u ) ) / _c;
}

Vector RadiationHydrodynamics1D::SF( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0]     = u[_nMoments + 1] / rho;
    v[1]     = 0.0;
    v[2]     = 0.0;
    double p = ( _gamma - 1.0 ) * ( u[_nMoments + 2] - 0.5 * rho * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    double T = p / ( _R * rho );
    return Fr0( u ) * ( -_sigmaT ) + v * _sigmaA * ( std::pow( T, 4 ) - Er ) / _c;
}

Matrix RadiationHydrodynamics1D::F( const Vector& u ) const {
    std::cerr << "F not tested" << std::endl;
    exit( EXIT_FAILURE );
    Matrix flux( u.size(), 2 );
    double rhoInv = 1.0 / u[_nMoments + 0];
    double v1     = u[_nMoments + 1] * rhoInv;
    double v2     = 0.0;
    double p      = ( _gamma - 1.0 ) * ( u[_nMoments + 2] - 0.5 * u[_nMoments + 0] * ( pow( v1, 2 ) + pow( v2, 2 ) ) );

    Vector momentFluxX = _Ax * u;
    Vector momentFluxY = _Ay * u;
    for( unsigned i = 0; i < _nMoments; ++i ) {
        flux( i, 0 ) = momentFluxX[i];
        flux( i, 1 ) = momentFluxY[i];
    }
    flux( _nMoments + 0, 0 ) = u[_nMoments + 1];
    flux( _nMoments + 1, 0 ) = u[_nMoments + 1] * v1 + p;
    flux( _nMoments + 2, 0 ) = 0.0;
    flux( _nMoments + 2, 0 ) = ( u[_nMoments + 2] + p ) * v1;
    flux( _nMoments + 0, 1 ) = u[_nMoments + 2];
    flux( _nMoments + 1, 1 ) = u[_nMoments + 2] * v1;
    flux( _nMoments + 2, 1 ) = u[_nMoments + 2] * v2 + p;
    flux( _nMoments + 2, 1 ) = ( u[_nMoments + 2] + p ) * v2;
    return flux;
}

Vector RadiationHydrodynamics1D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    /*
    double x0      = 0.0;
    double y0      = 0.0;
    double s2      = 3.2 * std::pow( 0.01, 2 );    // std::pow( 0.03, 2 );
    double s2Euler = 3.2 * std::pow( 0.15, 2 );    // std::pow( 0.03, 2 );
    double floor   = 0.0;
    _sigma         = _settings->GetSigma();

    y[0] = Delta( 0, 0 ) *
           std::fmax( floor, 1.0 / ( 4.0 * M_PI * s2 ) * exp( -( ( x[0] - x0 ) * ( x[0] - x0 ) + ( x[1] - y0 ) * ( x[1] - y0 ) ) / 4.0 / s2 ) );
    y[_nMoments + 0] = 1.0;
    y[_nMoments + 0] = std::fmax(
        floor, 1.0 / ( 4.0 * M_PI * s2Euler ) * exp( -( ( x[0] - x0 ) * ( x[0] - x0 ) + ( x[1] - y0 ) * ( x[1] - y0 ) ) / 4.0 / s2Euler ) );
    y[_nMoments+2] = std::fmax(
        floor, 1.0 / ( 4.0 * M_PI * s2Euler ) * exp( -( ( x[0] - x0 ) * ( x[0] - x0 ) + ( x[1] - y0 ) * ( x[1] - y0 ) ) / 4.0 / s2Euler ) );
    */
    if( x[0] < 0.5 ) {
        y[0]             = 1.0;
        y[_nMoments + 0] = 1.0;
        y[_nMoments + 2] = 1.0;
    }
    else {
        y[0]             = 1.0;
        y[_nMoments + 0] = 1.0;
        y[_nMoments + 2] = 1.0;
    }
    return y;
}
/*
void RadiationHydrodynamics1D::DS( Vector& ds, const Vector& u ) const {
    double gamma      = _gamma;
    double rho        = u[_nMoments + 0];
    double rhoU       = u[_nMoments + 1];
    double rhoV       = u[_nMoments + 2];
    double rhoV2      = pow( rhoV, 2 );
    double rhoU2      = pow( rhoU, 2 );
    double rhoE       = u[_nMoments+2];
    ds[_nMoments + 0] = ( rhoU2 + rhoV2 + gamma * ( 2 * rho * rhoE - rhoU2 - rhoV2 ) ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) -
                        std::log( pow( rho, gamma ) * ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );
    ds[_nMoments + 1] = -( ( 2 * rho * rhoU ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[_nMoments + 2] = -( ( 2 * rho * rhoV ) / ( -2 * rho * rhoE + rhoU2 + rhoV2 ) );
    ds[_nMoments+2] = -( rho / ( rhoE - ( rhoU2 + rhoV2 ) / ( 2 * rho ) ) );

    // safety factor on ds3
    // ds[3] -= 1e-7;
}*/

Matrix RadiationHydrodynamics1D::FRadiation( const Vector& u ) const {
    Matrix flux( u.size(), 1 );

    column( flux, 0 ) = _Ax * u;
    // column( flux, 1 ) = _Ay * u;
    // column( flux, 0 ) = _Az * u;

    return flux;
}

Matrix RadiationHydrodynamics1D::FEuler( const Vector& u ) const {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v + p;
    flux( 2, 0 ) = ( u[2] + p ) * v;

    return flux;
}

Vector RadiationHydrodynamics1D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    Vector out( _nStates, 0.0 );
    Vector outEuler( 3, false );
    Vector uRadiation( _nStates, false );
    Vector vRadiation( _nStates, false );
    Vector uEuler( _nStates, false );
    Vector vEuler( _nStates, false );
    // copy radiation part on uRadiation
    for( unsigned i = 0; i < _nMoments; ++i ) {
        uRadiation[i] = u[i];
        vRadiation[i] = v[i];
    }
    for( unsigned i = 0; i < 3; ++i ) {
        uEuler[i] = u[_nMoments + i];
        vEuler[i] = v[_nMoments + i];
    }
    double rhoInv = 1.0 / u[_nMoments + 0];
    double vU     = u[_nMoments + 1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[_nMoments + 2] - 0.5 * u[_nMoments + 0] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[_nMoments + 0];
    double vV     = v[_nMoments + 1] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[_nMoments + 2] - 0.5 * v[_nMoments + 0] * pow( vV, 2 ) );
    double aV     = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * vU;
    double uVProjected = nUnit[0] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        outEuler = FEuler( uEuler ) * n;
    else if( lambdaMax <= 0 )
        outEuler = FEuler( vEuler ) * n;
    else {
        outEuler = ( 1.0 / ( lambdaMax - lambdaMin ) ) *
                   ( lambdaMax * FEuler( uEuler ) * n - lambdaMin * FEuler( vEuler ) * n + lambdaMax * lambdaMin * ( vEuler - uEuler ) * norm( n ) );
    }

    // write radiation part on _nMoments entries
    out = FRadiation( 0.5 * ( uRadiation + vRadiation ) ) * n - 0.5 * ( vRadiation - uRadiation ) * norm( n );

    // save Euler part on return vector
    for( unsigned s = 0; s < 3; ++s ) out[_nMoments + s] = outEuler[s];    // TODO: DEBUG

    return out;
}

Matrix RadiationHydrodynamics1D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

double RadiationHydrodynamics1D::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    double dtMinTotal = 1e10;
    double dtMin;
    double rhoInv, uU, vU, p, a, cfl;
    unsigned kEnd = _settings->GetNqPEAtRef( level );

    cfl = _settings->GetCFL();

    for( unsigned k = 0; k < kEnd; ++k ) {
        rhoInv = 1.0 / u( _nMoments + 0, k );
        uU     = u( _nMoments + 1, k ) * rhoInv;
        vU     = 0.0;
        p      = ( _gamma - 1.0 ) * ( u( _nMoments + 2, k ) - 0.5 * u( _nMoments + 0, k ) * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
        a      = sqrt( _gamma * p * rhoInv );

        dtMin      = ( cfl / dx ) * std::min( std::min( std::fabs( 1.0 / ( vU - a ) ), std::fabs( 1.0 / ( vU + a ) ) ),
                                         std::min( std::fabs( 1.0 / ( uU + a ) ), std::fabs( 1.0 / ( uU - a ) ) ) );
        dtMinTotal = std::min( dtMin, dtMinTotal );
    }
    // std::cout << "P_N dt = " << _c * cfl / dx << ", Euler dt =  " << dtMinTotal << std::endl;
    return std::min( dtMinTotal, _c * cfl / dx );
}

Matrix RadiationHydrodynamics1D::BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const {
    unsigned nStates = u.rows();
    Vector outEuler( 3, false );
    Matrix flux( nStates, 2 );
    unsigned Nq = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    Vector uB( nStates );
    for( unsigned k = 0; k < Nq; ++k ) {
        // part radiative transfer
        Vector uM          = column( u, k );
        Vector correctTerm = nUnit[0] * _Ax * uM + nUnit[1] * _Ay * uM;
        column( flux, 0 )  = _Ax * uM - correctTerm * nUnit[0];
        column( flux, 1 )  = _Ay * uM - correctTerm * nUnit[1];
        column( y, k )     = flux * n;

        // part Euler
        Vector v( 2, 0.0 );
        v.reset();
        v[0]           = u( _nMoments + 1, k ) / u( _nMoments + 0, k );
        v[1]           = 0.0;
        double vn      = dot( nUnit, v );
        Vector Vn      = vn * nUnit;
        Vector Vb      = -Vn + v;
        double velMagB = Vb[0] * Vb[0] + Vb[1] * Vb[1];
        double velMag  = v[0] * v[0] + v[1] * v[1];
        double rho     = u( _nMoments + 0, k );
        uB[0]          = rho;
        uB[1]          = rho * ( Vb[0] );
        uB[2]          = rho * ( Vb[1] );
        uB[3]          = u( _nMoments + 2, k ) + rho * 0.5 * ( velMagB - velMag );
        outEuler       = FEuler( uB ) * n;

        // save Euler part on return vector
        for( unsigned s = 0; s < 3; ++s ) y( _nMoments + s, k ) = outEuler[s];    // TODO: DEBUG
    }
    return y;
}
