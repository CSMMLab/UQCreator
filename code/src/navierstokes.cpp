#include "navierstokes.h"

NavierStokes::NavierStokes( Settings* settings ) : Problem( settings ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Euler] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

NavierStokes::~NavierStokes() {}

Vector NavierStokes::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n, const double& length ) {

    // gas property
    double inK = (3.0 - _gamma) / (_gamma - 1.0);
    double vis = 0.001;
    
    // left interface
    double wL[3], primL[3];
    
    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];
    
    get_primitive(primL, wL, _gamma);
    
    double pL = ( _gamma - 1.0 ) * ( wL[2] - 0.5 * wL[0] * pow( wL[1], 2 ) );
    double ssL = sqrt( _gamma * pL / primL[0] );
    
    // right interface
    double wR[3], primR[3];
    
    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];
    
    get_primitive(primR, wR, _gamma);
    
    double pR = ( _gamma - 1.0 ) * ( wR[2] - 0.5 * wR[0] * pow( wR[1], 2 ) );
    double ssR = sqrt( _gamma * pR / primR[0] );

    // projection
    double uUProjected = nUnit[0] * primL[1];
    double uVProjected = nUnit[0] * primR[1];

    double lambdaMin = uUProjected - ssL;
    double lambdaMax = uVProjected + ssR;

    if( lambdaMin >= 0 )
        return F( u ) * nUnit;
    else if( lambdaMax <= 0 )
        return F( v ) * nUnit;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * F( u ) * nUnit - lambdaMin * F( v ) * nUnit + lambdaMax * lambdaMin * ( v - u ) );
    }





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


