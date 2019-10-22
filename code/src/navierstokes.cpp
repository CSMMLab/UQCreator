#include "navierstokes.h"
//#include "gks.cpp"

/* Transform between conservative and primitive variables */
void get_conserved( double* w, double prim[3], double gam ) {
    double smv = 1e-6;

    w[0] = prim[0];
    w[1] = prim[0] * prim[1];
    w[2] = 0.5 * prim[0] / ( prim[2] + smv ) / ( gam - 1.0 ) + 0.5 * prim[0] * prim[1] * prim[1];
}

void get_primitive( double* prim, double w[3], double gam ) {
    double smv = 1e-6;

    prim[0] = w[0];
    prim[1] = w[1] / ( w[0] + smv );
    prim[2] = 0.5 * w[0] / ( gam - 1.0 ) / ( w[2] - 0.5 * w[1] * w[1] / ( w[0] + smv ) + smv );
}

/* Calculate particle collision time */
double
get_tau( double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double mu, double dt ) {
    double smv = 1e-6;

    double tau_num = 10.0 * fabs( density_left / ( lambda_left + smv ) - density_right / ( lambda_right + smv ) ) /
                     fabs( density_left / ( lambda_left + smv ) + density_right / ( lambda_right + smv ) ) * dt;
    double tau_ns = 2.0 * mu * lambda0 / ( density0 + smv );

    return tau_num + tau_ns;
}

/* Calculate moments over particle velocity space */
void calc_moment( double* Mu, double* Mu_L, double* Mu_R, double* Mxi, double prim[3], double inK ) {

    double smv      = 1e-6;
    const double PI = 3.1415926535;

    // moments of normal velocity
    Mu_L[0] = 0.5 * erfc( -sqrt( prim[2] ) * prim[1] );
    Mu_L[1] = prim[1] * Mu_L[0] + 0.5 * exp( -prim[2] * pow( prim[1], 2 ) ) / sqrt( PI * prim[2] + smv );
    Mu_R[0] = 0.5 * erfc( sqrt( prim[2] ) * prim[1] );
    Mu_R[1] = prim[1] * Mu_R[0] - 0.5 * exp( -prim[2] * pow( prim[1], 2 ) ) / sqrt( PI * prim[2] + smv );

    for( int i = 2; i <= 6; i++ ) {
        Mu_L[i] = prim[1] * Mu_L[i - 1] + 0.5 * ( i - 1 ) * Mu_L[i - 2] / ( prim[2] + smv );
        Mu_R[i] = prim[1] * Mu_R[i - 1] + 0.5 * ( i - 1 ) * Mu_R[i - 2] / ( prim[2] + smv );
    }

    for( int i = 0; i <= 6; i++ ) {
        Mu[i] = Mu_L[i] + Mu_R[i];
    }

    // moments of inner degrees of freedom
    Mxi[0] = 1.0;
    Mxi[1] = 0.5 * inK / ( prim[2] + smv );
    Mxi[2] = ( pow( inK, 2 ) + 2.0 * inK ) / ( 4.0 * pow( prim[2], 2 ) + smv );
}

void moment_uv( double* moment_uv, double Mu[7], double Mxi[3], int alpha, int delta ) {
    moment_uv[0] = Mu[alpha] * Mxi[delta / 2];
    moment_uv[1] = Mu[alpha + 1] * Mxi[delta / 2];
    moment_uv[2] = 0.5 * ( Mu[alpha + 2] * Mxi[delta / 2] + Mu[alpha] * Mxi[( delta + 2 ) / 2] );
}

void moment_au( double* moment_au, double a[3], double Mu[7], double Mxi[3], int alpha ) {
    double t0[3], t1[3], t2[3], t3[3];

    moment_uv( t0, Mu, Mxi, alpha + 0, 0 );
    moment_uv( t1, Mu, Mxi, alpha + 1, 0 );
    moment_uv( t2, Mu, Mxi, alpha + 2, 0 );
    moment_uv( t3, Mu, Mxi, alpha + 0, 2 );

    for( int i = 0; i <= 2; i++ ) {
        moment_au[i] = a[0] * t0[i] + a[1] * t1[i] + 0.5 * a[2] * t2[i] + 0.5 * a[2] * t3[i];
    }
}

/* Gas kinetic flux function */
void GKS( double* flux, const Vector& u, const Vector& v, double gam, double mu, double dt ) {

    /* HLL Flux */
    // left interface
    double wL[3], primL[3];

    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];

    get_primitive( primL, wL, gam );

    double pL  = ( gam - 1.0 ) * ( wL[2] - 0.5 * wL[0] * pow( wL[1], 2 ) );
    double ssL = sqrt( gam * pL / wL[0] );

    // right interface
    double wR[3], primR[3];

    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];

    get_primitive( primR, wR, gam );

    double pR  = ( gam - 1.0 ) * ( wR[2] - 0.5 * wR[0] * pow( wR[1], 2 ) );
    double ssR = sqrt( gam * pR / wR[0] );

    double lambdaMin = primL[1] - ssL;
    double lambdaMax = primR[1] + ssR;

    if( lambdaMin >= 0 ) {
        flux[0] = wL[1];
        flux[1] = wL[1] * primL[1] + pL;
        flux[2] = ( wL[2] + pL ) * primL[1];
    }
    else if( lambdaMax <= 0 ) {
        flux[0] = wR[1];
        flux[1] = wR[1] * primR[1] + pR;
        flux[2] = ( wR[2] + pR ) * primR[1];
    }
    else {
        flux[0] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * wL[1] - lambdaMin * wR[1] + lambdaMax * lambdaMin * ( wR[0] - wL[0] ) );
        flux[1] = ( 1.0 / ( lambdaMax - lambdaMin ) ) *
                  ( lambdaMax * ( wL[1] * primL[1] + pL ) - lambdaMin * ( wR[1] * primR[1] + pR ) + lambdaMax * lambdaMin * ( wR[1] - wL[1] ) );
        flux[2] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * ( ( wL[2] + pL ) * primL[1] ) - lambdaMin * ( ( wR[2] + pR ) * primR[1] ) +
                                                          lambdaMax * lambdaMin * ( wR[2] - wL[2] ) );
    }

    /* GKS
    // gas property
    double inK = ( 3.0 - gam ) / ( gam - 1.0 );

    // left interface
    double wL[3], primL[3];

    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];

    get_primitive( primL, wL, gam );

    // right interface
    double wR[3], primR[3];

    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];

    get_primitive( primR, wR, gam );

    // central interface
    double Mu[7], MuL[7], MuR[7], Mxi[3];
    double Mau[3], MauL[3], MauR[3];

    calc_moment( Mu, MuL, MuR, Mxi, primL, inK );
    moment_uv( MauL, MuL, Mxi, 0, 0 );

    calc_moment( Mu, MuL, MuR, Mxi, primR, inK );
    moment_uv( MauR, MuR, Mxi, 0, 0 );

    double w[3], prim[3], tau;
    for( int i = 0; i <= 2; i++ ) {
        w[i] = primL[0] * MauL[i] + primR[0] * MauR[i];
    }

    get_primitive( prim, w, gam );

    tau = get_tau( primL[0], primR[0], prim[0], primL[2], primR[2], prim[2], mu, dt );

    // time integral terms
    double Mt[5];

    Mt[3] = dt;//tau * ( 1.0 - exp( -dt / tau ) );
    //Mt[4] = -tau * dt * exp( -dt / tau ) + tau * Mt[3];
    Mt[0] = dt - Mt[3];
    //Mt[1] = -tau * Mt[0] + Mt[4];
    //Mt[2] = dt * dt / 2.0 - tau * Mt[0];

    // calculate the flux of conservative variables related to g0
    calc_moment( Mu, MuL, MuR, Mxi, prim, inK );
    moment_uv( Mau, Mu, Mxi, 1, 0 );

    for( int i = 0; i <= 2; i++ ) {
        flux[i] = 0;//Mt[0] * prim[0] * Mau[i];
    }

    // calculate the flux of conservative variables related to f0
    calc_moment( Mu, MuL, MuR, Mxi, primL, inK );
    moment_uv( MauL, MuL, Mxi, 1, 0 );
    for( int i = 0; i <= 2; i++ ) {
        flux[i] += Mt[3] * primL[0] * MauL[i];
    }

    calc_moment( Mu, MuL, MuR, Mxi, primR, inK );
    moment_uv( MauR, MuR, Mxi, 1, 0 );
    for( int i = 0; i <= 2; i++ ) {
        flux[i] += Mt[3] * primR[0] * MauR[i];
    }

    // final flux
    for( int i = 0; i <= 2; i++ ) {
        flux[i] = flux[i] / dt;
    }
    */
}

NavierStokes::NavierStokes( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
        _mu = problem->get_as<double>( "mu" ).value_or( 0.01 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Euler] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

NavierStokes::~NavierStokes() {}

Vector NavierStokes::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    /* HLL Flux */
    /*
    // left interface
    double wL[3], primL[3];
    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];
    get_primitive( primL, wL, _gamma );
    double pL  = ( _gamma - 1.0 ) * ( wL[2] - 0.5 * primL[0] * pow( primL[1], 2 ) );
    double ssL = sqrt( _gamma * pL / wL[0] );

    // right interface
    double wR[3], primR[3];
    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];
    get_primitive( primR, wR, _gamma );
    double pR  = ( _gamma - 1.0 ) * ( wR[2] - 0.5 * primR[0] * pow( primR[1], 2 ) );
    double ssR = sqrt( _gamma * pR / wR[0] );

    // projection
    double uUProjected = nUnit[0] * primL[1];
    double uVProjected = nUnit[0] * primR[1];
    double lambdaMin = uUProjected - ssL;
    double lambdaMax = uVProjected + ssR;
    //double lambdaMin = primL[1] - ssL;
    //double lambdaMax = primR[1] + ssR;

    // central
    //return 0.5 * ( F(u) + F(v) ) * nUnit;

    // original
    //if( lambdaMin >= 0 )
    //    return F( u ) * nUnit;
    //else if( lambdaMax <= 0 )
    //    return F( v ) * nUnit;
    //else {
    //    return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * F( u ) * nUnit - lambdaMin * F( v ) * nUnit + lambdaMax * lambdaMin * ( v - u )
    );
    //}

    // modified 1
    ////std::cout<<"nUnit = "<<nUnit<<std::endl;
    //double interfaceSign = nUnit[0];

    //if( lambdaMin >= 0 )
    //    return FF( u ) * interfaceSign;
    //else if( lambdaMax <= 0 )
    //    return FF( v ) * interfaceSign;
    //else {
    //    return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * FF( u ) * interfaceSign - lambdaMin * FF( v ) * interfaceSign + lambdaMax *
    lambdaMin * ( v - u ) );
    //}

    // modified 2
    double interfaceSign = nUnit[0];
    Vector flux(3);
    if( lambdaMin > 0 ) {
        flux[0] = wL[1] * interfaceSign;
        flux[1] = ( wL[1] * primL[1] + pL ) * interfaceSign;
        flux[2] = ( (wL[2] + pL) * primL[1] ) * interfaceSign;
    }
    else if( lambdaMax < 0 ) {
        flux[0] = wR[1] * interfaceSign;
        flux[1] = ( wR[1] * primR[1] + pR ) * interfaceSign;
        flux[2] = ( (wR[2] + pR) * primR[1] ) * interfaceSign;
    }
    else {
        flux[0] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * wL[1] * interfaceSign - lambdaMin * wR[1] * interfaceSign + lambdaMax *
    lambdaMin * ( v[0] - u[0] ) ); flux[1] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * ( wL[1] * primL[1] + pL ) * interfaceSign - lambdaMin
    * ( wR[1] * primR[1] + pR ) * interfaceSign + lambdaMax * lambdaMin * ( v[1] - u[1] ) ); flux[2] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * (
    lambdaMax * ( (wL[2] + pL) * primL[1] ) * interfaceSign - lambdaMin * ( (wR[2] + pR) * primR[1] ) * interfaceSign + lambdaMax * lambdaMin * ( v[2]
    - u[2] ) );
    }

    return flux;
    */

    /* GKS Flux */

    // computational parameters
    double dx            = norm( n );
    double dt            = _settings->GetDT();
    double inK           = ( 3.0 - _gamma ) / ( _gamma - 1.0 );
    double interfaceSign = nUnit[0];
    double Mu[7], MuL[7], MuR[7], Mxi[3];
    double Mau[3], MauL[3], MauR[3];

    double pL, ssL, pR, ssR;
    double wR[3], primR[3], wL[3], primL[3];
    if( interfaceSign > 1 ) {
        // left interface
        wL[0] = u[0];
        wL[1] = u[1];
        wL[2] = u[2];
        get_primitive( primL, wL, _gamma );
        pL  = ( _gamma - 1.0 ) * ( wL[2] - 0.5 * primL[0] * pow( primL[1], 2 ) );
        ssL = sqrt( _gamma * pL / wL[0] );

        // right interface

        wR[0] = v[0];
        wR[1] = v[1];
        wR[2] = v[2];
        get_primitive( primR, wR, _gamma );
        pR  = ( _gamma - 1.0 ) * ( wR[2] - 0.5 * primR[0] * pow( primR[1], 2 ) );
        ssR = sqrt( _gamma * pR / wR[0] );
    }
    else {
        // left interface
        wL[0] = v[0];
        wL[1] = v[1];
        wL[2] = v[2];
        get_primitive( primL, wL, _gamma );
        pL  = ( _gamma - 1.0 ) * ( wL[2] - 0.5 * primL[0] * pow( primL[1], 2 ) );
        ssL = sqrt( _gamma * pL / wL[0] );

        // right interface

        wR[0] = u[0];
        wR[1] = u[1];
        wR[2] = u[2];
        get_primitive( primR, wR, _gamma );
        pR  = ( _gamma - 1.0 ) * ( wR[2] - 0.5 * primR[0] * pow( primR[1], 2 ) );
        ssR = sqrt( _gamma * pR / wR[0] );
    }

    // central interface
    // calc_moment( Mu, MuL, MuR, Mxi, primL, inK ); moment_uv( MauL, MuL, Mxi, 0, 0 );
    // calc_moment( Mu, MuL, MuR, Mxi, primR, inK ); moment_uv( MauR, MuR, Mxi, 0, 0 );

    // double w[3], prim[3], tau;
    // for( int i = 0; i <= 2; i++ ) {
    //    w[i] = primL[0] * MauL[i] + primR[0] * MauR[i];
    //}
    // get_primitive( prim, w, _gamma );
    // tau = get_tau( primL[0], primR[0], prim[0], primL[2], primR[2], prim[2], _mu, dt );

    // time integral terms
    double Mt[5];
    Mt[3] = dt;    // tau * ( 1.0 - exp( -dt / tau ) );
    // Mt[4] = -tau * dt * exp( -dt / tau ) + tau * Mt[3];
    // Mt[0] = 0.0;//dt - Mt[3];
    // Mt[1] = -tau * Mt[0] + Mt[4];
    // Mt[2] = dt * dt / 2.0 - tau * Mt[0];

    // calculate the flux of conservative variables related to g0
    // calc_moment( Mu, MuL, MuR, Mxi, prim, inK ); moment_uv( Mau, Mu, Mxi, 1, 0 );

    Vector flux( 3 );
    for( int i = 0; i < 3; i++ ) {
        flux[i] = 0.0;    // Mt[0] * prim[0] * Mau[i];
    }

    // calculate the flux of conservative variables related to f0
    calc_moment( Mu, MuL, MuR, Mxi, primL, inK );
    moment_uv( MauL, MuL, Mxi, 1, 0 );
    for( int i = 0; i < 3; i++ ) {
        flux[i] += Mt[3] * primL[0] * MauL[i];
    }

    calc_moment( Mu, MuL, MuR, Mxi, primR, inK );
    moment_uv( MauR, MuR, Mxi, 1, 0 );
    for( int i = 0; i < 3; i++ ) {
        flux[i] += Mt[3] * primR[0] * MauR[i];
    }

    // final flux
    for( int i = 0; i < 3; i++ ) {
        flux[i] = flux[i] / dt * interfaceSign;
    }

    return flux;
}

Matrix NavierStokes::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix NavierStokes::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v + p;
    flux( 2, 0 ) = ( u[2] + p ) * v;
    return flux;
}

Vector NavierStokes::FF( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Vector flux( u.size() );
    flux[0] = u[1];
    flux[1] = u[1] * v + p;
    flux[2] = ( u[2] + p ) * v;
    return flux;
}

Vector NavierStokes::IC( const Vector& x, const Vector& xi ) {
    double x0    = 0.5;
    double gamma = 1.4;

    double rhoL = 1.0;
    double rhoR = 0.9;
    double pL   = 1.0;
    double pR   = 0.9;
    double uL   = 0.0;
    double uR   = 0.0;
    Vector y( _nStates );
    _sigma = _settings->GetSigma();
    if( x[0] < x0 + _sigma[0] * xi[0] ) {
        y[0]                  = rhoL;
        y[1]                  = rhoL * uL;
        double kineticEnergyL = 0.5 * rhoL * pow( uL, 2 );
        double innerEnergyL   = ( pL / ( rhoL * ( gamma - 1 ) ) ) * rhoL;
        y[2]                  = kineticEnergyL + innerEnergyL;
    }
    else {
        y[0] = rhoR;
        if( xi.size() > 1 ) {
            y[0] += _sigma[1] * xi[1];
        }
        if( xi.size() > 2 ) {
            pR += _sigma[2] * xi[2];
        }
        y[1]                  = rhoR * uR;
        double kineticEnergyR = 0.5 * rhoR * pow( uR, 2 );
        double innerEnergyR   = ( pR / ( rhoR * ( gamma - 1 ) ) ) * rhoR;
        y[2]                  = kineticEnergyR + innerEnergyR;
    }
    return y;
}
