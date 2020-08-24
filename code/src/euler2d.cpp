#include "euler2d.h"

Euler2D::Euler2D( Settings* settings ) : Problem( settings ), _problemType( I_NACA ) {
    _nStates = 4;
    _settings->SetNStates( _nStates );
    _sigma = _settings->GetSigma();

    try {
        auto file     = cpptoml::parse_file( _settings->GetInputFile() );
        auto general  = file->get_table( "general" );
        auto ICString = general->get_as<std::string>( "testCase" );
        if( ICString ) {
            if( ICString->compare( "nozzle" ) == 0 ) {
                _problemType = ICType::I_NOZZLE;
            }
            else if( ICString->compare( "nozzleSod" ) == 0 ) {
                _problemType = ICType::I_NOZZLE_SOD;
            }
            else if( ICString->compare( "naca" ) == 0 ) {
                _problemType = ICType::I_NACA;
            }
            else if( ICString->compare( "nacaHighMach" ) == 0 ) {
                _problemType = ICType::I_NACA_HIGHMACH;
            }
            else {
                _log->error( "[euler2d] Unknown testcase defined!" );
            }
        }

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[euler2d] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Euler2D::~Euler2D() {}

void Euler2D::Solve() {}

Vector Euler2D::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    double rhoInv = 1.0 / u[0];
    double uU     = u[1] * rhoInv;
    double vU     = u[2] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
    double aU     = sqrt( _gamma * p * rhoInv );

    rhoInv    = 1.0 / v[0];
    double uV = v[1] * rhoInv;
    double vV = v[2] * rhoInv;
    p         = ( _gamma - 1.0 ) * ( v[3] - 0.5 * v[0] * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
    double aV = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
    double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        return F( u ) * n;
    else if( lambdaMax <= 0 )
        return F( v ) * n;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
}

Matrix Euler2D::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        // column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );

        double rhoInv = 1.0 / u( 0, k );
        double uU     = u( 1, k ) * rhoInv;
        double vU     = u( 2, k ) * rhoInv;
        double p      = ( _gamma - 1.0 ) * ( u( 3, k ) - 0.5 * u( 0, k ) * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
        double aU     = sqrt( _gamma * p * rhoInv );

        rhoInv    = 1.0 / v( 0, k );
        double uV = v( 1, k ) * rhoInv;
        double vV = v( 2, k ) * rhoInv;
        p         = ( _gamma - 1.0 ) * ( v( 3, k ) - 0.5 * v( 0, k ) * ( pow( uV, 2 ) + pow( vV, 2 ) ) );
        double aV = sqrt( _gamma * p * rhoInv );

        double uUProjected = nUnit[0] * uU + nUnit[1] * vU;
        double uVProjected = nUnit[0] * uV + nUnit[1] * vV;

        double lambdaMin = uUProjected - aU;
        double lambdaMax = uVProjected + aV;

        if( lambdaMin >= 0 )
            column( y, k ) = F( column( u, k ) ) * n;
        else if( lambdaMax <= 0 )
            column( y, k ) = F( column( v, k ) ) * n;
        else {
            column( y, k ) = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * F( column( u, k ) ) - lambdaMin * F( column( v, k ) ) ) * n +
                                                                     lambdaMax * lambdaMin * ( column( v, k ) - column( u, k ) ) * norm( n ) );
        }
    }
    return y;
}

Matrix Euler2D::F( const Vector& u ) const {
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

Matrix Euler2D::F( const Matrix& u ) {
    _log->error( "[euler2d] Flux not implemented" );
    exit( EXIT_FAILURE );
    return 0.5 * pow( u, 2 );
}

Matrix Euler2D::BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    Vector uB( nStates );
    for( unsigned k = 0; k < Nq; ++k ) {
        Vector v( 2, 0.0 );
        v.reset();
        v[0]           = u( 1, k ) / u( 0, k );
        v[1]           = u( 2, k ) / u( 0, k );
        double vn      = dot( nUnit, v );
        Vector Vn      = vn * nUnit;
        Vector Vb      = v - Vn;
        double velMagB = Vb[0] * Vb[0] + Vb[1] * Vb[1];
        double velMag  = v[0] * v[0] + v[1] * v[1];
        double rho     = u( 0, k );
        uB[0]          = rho;
        uB[1]          = rho * ( Vb[0] );
        uB[2]          = rho * ( Vb[1] );
        uB[3]          = u( 3, k ) + rho * 0.5 * ( velMagB - velMag );
        column( y, k ) = F( uB ) * n;
    }
    return y;
}

double Euler2D::ComputeDt( const Tensor& u, double dx, unsigned level ) const {
    double dtMinTotal = 1e10;
    double dtMin;
    double rhoInv, uU, vU, p, a, cfl;
    unsigned kEnd = _settings->GetNqPEAtRef( level );

    cfl = _settings->GetCFL();
    for( unsigned l = 0; l < _settings->GetNMultiElements(); ++l ) {
        for( unsigned k = 0; k < kEnd; ++k ) {
            rhoInv = 1.0 / u( 0, l, k );
            uU     = u( 1, l, k ) * rhoInv;
            vU     = u( 2, l, k ) * rhoInv;
            p      = ( _gamma - 1.0 ) * ( u( 3, l, k ) - 0.5 * u( 0, l, k ) * ( pow( uU, 2 ) + pow( vU, 2 ) ) );
            a      = sqrt( _gamma * p * rhoInv );

            dtMin      = ( cfl / dx ) * std::min( std::min( std::fabs( 1.0 / ( vU - a ) ), std::fabs( 1.0 / ( vU + a ) ) ),
                                             std::min( std::fabs( 1.0 / ( uU + a ) ), std::fabs( 1.0 / ( uU - a ) ) ) );
            dtMinTotal = std::min( dtMin, dtMinTotal );
        }
    }

    return dtMinTotal;
}

Vector Euler2D::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );

    if( _problemType == I_NOZZLE || _problemType == I_NOZZLE_SOD ) {
        double gamma = 1.4;

        double uF = 0.0;
        double vF = 0.0;

        double rhoL = 1.0;
        double pL   = 1.0;

        y[0] = rhoL;
        y[1] = rhoL * uF;
        y[2] = rhoL * vF;
        y[3] = pL / ( gamma - 1 );

        if( x[0] > -0.5 + _sigma[0] * xi[0] ) {
            double rhoR = 0.8;
            double pR   = 0.125;    // 0.3;
            if( _problemType == I_NOZZLE_SOD ) {
                rhoR = 0.125;
                pR   = 0.1;
            }
            y[0] = rhoR;
            y[1] = rhoR * uF;
            y[2] = rhoR * vF;
            y[3] = pR / ( gamma - 1 );
        }
        return y;
    }
    if( _problemType == I_NACA || _problemType == I_NACA_HIGHMACH ) {
        double gamma       = 1.4;
        double R           = 287.87;
        double T           = 273.15;
        double Ma          = 0.8;
        double AoA         = 1.25;
        double AoAScaling  = 1.0;
        double p           = 101325.0;
        double rhoFarfield = p / ( R * T );

        if( _problemType == I_NACA_HIGHMACH ) {
            Ma          = 20.0;    // 6;
            rhoFarfield = 0.001027;
            AoA         = 10.0;    // 5
            AoAScaling  = 0.0;
            if( xi.size() == 1 ) {
                rhoFarfield = rhoFarfield + _sigma[0];
                rhoFarfield = rhoFarfield + xi[0] * _sigma[0];
            }
            p = rhoFarfield * ( R * T );
        }

        if( xi.size() > 1 ) {
            Ma = Ma - _sigma[1];
            Ma = Ma + xi[1] * _sigma[1];
        }

        double a = sqrt( gamma * R * T );

        double uMax  = Ma * a;
        double angle = ( AoA + AoAScaling * _sigma[0] * xi[0] ) * ( 2.0 * M_PI ) / 360.0;
        double uF    = uMax * cos( angle );
        double vF    = uMax * sin( angle );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        return y;
    }
}

Vector Euler2D::LoadIC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );

    double rhoFarfield = x[0];
    double u           = x[1] / rhoFarfield;
    double v           = x[2] / rhoFarfield;
    double uMax        = std::sqrt( u * u + v * v );
    double angle       = ( 1.25 + _sigma[0] * xi[0] ) * ( 2.0 * M_PI ) / 360.0;
    double uF          = uMax * cos( angle );
    double vF          = uMax * sin( angle );

    double gamma = 1.4;
    double p     = 101325.0;

    y[0]                  = rhoFarfield;
    y[1]                  = rhoFarfield * uF;
    y[2]                  = rhoFarfield * vF;
    double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
    double innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
    y[3]                  = kineticEnergyL + innerEnergyL;
    return y;
}
