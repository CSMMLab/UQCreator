#include "navierstokes.h"
#include "gks.cpp"

NavierStokes::NavierStokes( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
        _vis = problem->get_as<double>( "vis" ).value_or( 0.001 );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Euler] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

NavierStokes::~NavierStokes() {}

Vector NavierStokes::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    // gas property
    double inK = ( 3.0 - _gamma ) / ( _gamma - 1.0 );
    double dx  = norm( n );
    //double vis = _vis;

    // left interface
    double wL[3], primL[3];

    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];

    get_primitive( primL, wL, _gamma );

    double pL  = ( _gamma - 1.0 ) * ( wL[2] - 0.5 * wL[0] * pow( wL[1], 2 ) );
    double ssL = sqrt( _gamma * pL / primL[0] );

    // right interface
    double wR[3], primR[3];

    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];

    get_primitive( primR, wR, _gamma );

    double pR  = ( _gamma - 1.0 ) * ( wR[2] - 0.5 * wR[0] * pow( wR[1], 2 ) );
    double ssR = sqrt( _gamma * pR / primR[0] );

    // projection
    //double uUProjected = nUnit[0] * primL[1];
    //double uVProjected = nUnit[0] * primR[1];

    //double lambdaMin = uUProjected - ssL;
    //double lambdaMax = uVProjected + ssR;

    double dt = _settings->GetDT();

    // central interface
    double Mu[7], MuL[7], MuR[7], Mxi[3];
    double Mau[3], MauL[3], MauR[3];

    calc_moment(Mu, MuL, MuR, Mxi, primL, inK);
    moment_uv(MauL, Mu, Mxi, 0, 0);
    
    calc_moment(Mu, MuL, MuR, Mxi, primR, inK);
    moment_uv(MauR, Mu, Mxi, 0, 0);

    double w[3], prim[3], tau;
    for (int i=0;i<=2;i++)
    {
        w[i] = primL[0] * MauL[i] + primR[1] * MauR[i];
    }

    get_primitive(prim, w, gam);
    
    tau = get_tau(primL[0], primR[0], prim[0], primL[2], primR[2], prim[2], mu, dt);
    
    // time integral terms
    double Mt[5];

    Mt[3] = tau * (1.0 - exp(-dt / tau));
    Mt[4] = -tau * dt * exp(-dt / tau) + tau * Mt[3];
    Mt[0] = dt - Mt[3];
    Mt[1] = -tau * Mt[0] + Mt[4];
    Mt[2] = dt * dt / 2.0 - tau * Mt[0];

    // calculate the flux of conservative variables related to g0
    calc_moment(Mu, MuL, MuR, Mxi, prim, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    
    double flux[3]
    for (int i=0;i<=2;i++)
    {
        flux[i] = Mt[0] * primL[0] * Mau[i];
    }

    // calculate the flux of conservative variables related to f0
    calc_moment(Mu, MuL, MuR, Mxi, primL, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    for (int i=0;i<=2;i++)
    {
        flux[i] += Mt[3] * primL[0] * Mau[i];
    }
    
    calc_moment(Mu, MuL, MuR, Mxi, primR, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    for (int i=0;i<=2;i++)
    {
        flux[i] += Mt[3] * primR[0] * Mau[i];
    }
    
    // final flux
    for (int i=0;i<=2;i++)
    {
        flux[i] = flux[i] / dt;
    }

    // return value
    Matrix fluxMatrix( u.size(), 1 );
    fluxMatrix( 0, 0 ) = flux[0];
    fluxMatrix( 1, 0 ) = flux[1];
    fluxMatrix( 2, 0 ) = flux[2];

    return fluxMatrix;

    //if( lambdaMin >= 0 )
    //    return F( u ) * nUnit;
    //else if( lambdaMax <= 0 )
    //    return F( v ) * nUnit;
    //else {
    //    return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * F( u ) * nUnit - lambdaMin * F( v ) * nUnit + lambdaMax * lambdaMin * ( v - u ) );
    //}
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

Vector NavierStokes::IC( const Vector& x, const Vector& xi ) {
    double x0    = 0.3;
    double gamma = 1.4;

    double rhoL = 1.0;
    double rhoR = 0.3;
    double pL   = 1.0;
    double pR   = 0.3;
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
