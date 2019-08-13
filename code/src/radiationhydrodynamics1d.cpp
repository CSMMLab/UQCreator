#include "radiationhydrodynamics1d.h"

RadiationHydrodynamics1D::RadiationHydrodynamics1D( Settings* settings ) : PNEquations( settings, true ) {
    _N = 1;    // set moment order

    // specify test case
    _testCase = 2;

    // constants
    double k = 1.38064852e-16;         // Boltzmann's constant [cm^2 g / (s^2 K)]
    _gamma   = 5.0 / 3.0;              // adiabatic constant
    _R       = 8.3144621 * 1e7;        // specific gas constant [erg / (K mol)]
    _cLight  = 299792458.0 * 100.0;    // speed of light in [cm/s]
    _aR      = 7.5657 * 1e-15;         // radiation constant [erg/(cm^3 K^4)]
    _h       = 6.62607004 * 1e-27;     // Planck's constant [cm^2 g/s]

    // initial chock states
    double rhoL, uL, TL, rhoR, uR, TR, lRef;
    if( _testCase == 1 ) {
        rhoL = 5.45887 * 1e-13;    // [g/cm^3]
        uL   = 2.3545 * 1e5;       // [cm/s]
        TL   = 100.0;              // [K]
        rhoR = 1.2479 * 1e-12;
        uR   = 1.03 * 1e5;
        TR   = 207.757;
        lRef = 1000;
    }
    else {
        rhoL = 5.45887 * 1e-13;    // [g/cm^3]
        uL   = 2.3545 * 1e5;       // [cm/s]
        TL   = 100.0;              // [K]
        rhoR = 1.964050 * 1e-12;
        uR   = 1.63 * 1e5;
        TR   = 855.72;
        lRef = 4000;
    }

    // reference values
    _rhoRef = rhoR;                           // reference density [g/cm^3]
    _TRef   = TR;                             // reference temperature [K]
    _aRef   = sqrt( _gamma * _R * _TRef );    // reference speed of sound [cm/s]
    _pRef   = _rhoRef * pow( _aRef, 2 );      // update pRef with speed of sound

    std::cout << "vL = " << uL << ", aRef = " << _aRef << std::endl;

    _c = _cLight / _aRef;
    _P = _aR * pow( _TRef, 4 ) / ( _rhoRef * pow( _aRef, 2 ) );    // * 100;

    // set cross sections
    _sigmaA = 3.93 * 1e-5 * lRef;    // absorption coefficient
    _sigmaT = 0.848902 * lRef;       // 7.5657 * 1e-15;
    _sigmaS = _sigmaT - _sigmaA;     // scattering coefficient

    std::cout << "C = " << _c << " , P = " << _P << std::endl;

    _settings->SetGamma( _gamma );

    if( _N == 1 )    // GlobalIndex has different ordering here
        _nMoments = 4;
    else
        _nMoments = unsigned( GlobalIndex( _N, _N ) + 1 );
    _nStates = _nMoments + 3;    // total number of states in equations in P_N + Euler
    _settings->SetNStates( _nStates );
    _settings->SetSource( true );    // TODO:DEBUGGING
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
        SetupSystemMatrices();
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[pnequations] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
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
                _Ax( i, n ) = _Ax( i, n ) * Delta( l, k );
                _Ay( i, n ) = _Ay( i, n ) * Delta( l, k );
                _Az( i, n ) = _Az( i, n ) * Delta( l, k );
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

double RadiationHydrodynamics1D::SE( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0]     = u[_nMoments + 1] / rho;
    v[1]     = 0.0;
    v[2]     = 0.0;
    double T = _gamma * ( _gamma - 1.0 ) * ( u[_nMoments + 2] / rho - 0.5 * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    return _sigmaA * ( std::pow( T, 4 ) - Er ) + ( _sigmaA - _sigmaS ) * v.inner( Fr0( u ) ) / _c;
}

Vector RadiationHydrodynamics1D::SF( const Vector& u ) const {
    double Er  = u[0];
    double rho = u[_nMoments + 0];
    Vector v( 3, false );
    v[0] = u[_nMoments + 1] / rho;
    v[1] = 0.0;
    v[2] = 0.0;
    // double p = _gamma * ( _gamma - 1.0 ) * ( u[_nMoments + 2] - 0.5 * rho * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );
    double T = _gamma * ( _gamma - 1.0 ) * ( u[_nMoments + 2] / rho - 0.5 * ( pow( v[0], 2 ) + pow( v[1], 2 ) + pow( v[2], 2 ) ) );

    return Fr0( u ) * ( -_sigmaT ) + v * _sigmaA * ( std::pow( T, 4 ) - Er ) / _c;
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
    Matrix Pr( _nMoments - 1, _nMoments - 1, false );
    Vector v0Pr = v[0] * _Ax * u;
    Vector v1Pr = v[1] * _Ay * u;
    Vector v2Pr = v[2] * _Az * u;
    Vector vPr( 3, false );

    // get rid of 0th entry
    for( unsigned s = 0; s < 3; ++s ) vPr[s] = v0Pr[s + 1] + v1Pr[s + 1] + v2Pr[s + 1];

    return Fr - ( v * Er + vPr ) / _c;
}

Vector RadiationHydrodynamics1D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double k = 1.38064852e-16;    // Boltzmann's constant [cm^2 g / (s^2 K)]

    // initial chock states
    double rhoL, uL, TL, rhoR, uR, TR;
    rhoL = 5.45887 * 1e-13;    // [g/cm^3]
    uL   = 2.3545 * 1e5;       // [cm/s]
    TL   = 100.0;              // [K]
    if( _testCase == 1 ) {
        rhoR = 1.2479 * 1e-12;
        uR   = 1.03 * 1e5;
        TR   = 207.757;
    }
    else {
        rhoR = 1.964050 * 1e-12;
        uR   = 1.63 * 1e5;
        TR   = 855.72;
    }

    double shockPosition = 0.0;
    shockPosition += _settings->GetSigma()[0] * xi[0];

    if( x[0] < shockPosition ) {
        y[_nMoments + 0]      = rhoL / _rhoRef;
        y[_nMoments + 1]      = rhoL * uL / ( _rhoRef * _aRef );
        double pL             = TL * ( _R * rhoL ) / _pRef;
        double kineticEnergyL = 0.5 * rhoL * pow( uL, 2 ) / ( _rhoRef * pow( _aRef, 2 ) );
        double innerEnergyL   = ( pL / ( _gamma - 1.0 ) );
        y[_nMoments + 2]      = kineticEnergyL + innerEnergyL;

        double rho = y[_nMoments + 0];
        double v   = y[_nMoments + 1] / rho;
        y[0]       = pow( _gamma * ( _gamma - 1.0 ) * ( y[_nMoments + 2] / rho - 0.5 * ( pow( v, 2 ) ) ), 4 );
    }
    else {
        y[_nMoments + 0]      = rhoR / _rhoRef;
        y[_nMoments + 1]      = rhoR * uR / ( _rhoRef * _aRef );
        double pR             = TR * ( _R * rhoR ) / _pRef;
        double kineticEnergyR = 0.5 * rhoR * pow( uR, 2 ) / ( _rhoRef * pow( _aRef, 2 ) );
        double innerEnergyR   = ( pR / ( _gamma - 1.0 ) );
        y[_nMoments + 2]      = kineticEnergyR + innerEnergyR;

        double rho = y[_nMoments + 0];
        double v   = y[_nMoments + 1] / rho;
        y[0]       = pow( _gamma * ( _gamma - 1.0 ) * ( y[_nMoments + 2] / rho - 0.5 * ( pow( v, 2 ) ) ), 4 );
    }
    return y;
}

Matrix RadiationHydrodynamics1D::FRadiation( const Vector& u ) const {
    Matrix flux( u.size(), 1 );

    column( flux, 0 ) = _c * _Ax * u;
    // column( flux, 1 ) = _c * _Ay * u;
    // column( flux, 0 ) = _c * _Az * u;

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
    double dx = norm( n );
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
        outEuler = FEuler( uEuler ) * nUnit;
    else if( lambdaMax <= 0 )
        outEuler = FEuler( vEuler ) * nUnit;
    else {
        outEuler = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * FEuler( uEuler ) * nUnit - lambdaMin * FEuler( vEuler ) * nUnit +
                                                           lambdaMax * lambdaMin * ( vEuler - uEuler ) * norm( nUnit ) );
    }

    // write radiation part on _nMoments entries
    // std::cout << "c = " << _c << ", LF = " << 1 / ( dx / _settings->GetDT() ) << std::endl;
    // out = FRadiation( 0.5 * ( uRadiation + vRadiation ) ) * nUnit - 0.5 * ( vRadiation - uRadiation ) * dx / _settings->GetDT();
    out = FRadiation( 0.5 * ( uRadiation + vRadiation ) ) * nUnit - 0.5 * _c * ( vRadiation - uRadiation ) * norm( nUnit );

    // save Euler part on return vector
    for( unsigned s = 0; s < 3; ++s ) out[_nMoments + s] = outEuler[s];

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

    // std::cout << "dx = " << dx << std::endl;

    for( unsigned k = 0; k < kEnd; ++k ) {
        rhoInv = 1.0 / u( _nMoments + 0, k );
        uU     = u( _nMoments + 1, k ) * rhoInv;
        vU     = 0.0;
        p      = ( _gamma - 1.0 ) * ( u( _nMoments + 2, k ) - 0.5 * u( _nMoments + 0, k ) * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
        a      = sqrt( _gamma * p * rhoInv );

        dtMin      = ( cfl * dx ) * std::min( std::fabs( 1.0 / ( uU - a ) ), std::fabs( 1.0 / ( uU + a ) ) );
        dtMinTotal = std::min( dtMin, dtMinTotal );
    }
    return std::min( dtMinTotal, cfl * dx / _c );
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
