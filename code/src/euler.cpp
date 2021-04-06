#include "euler.h"

Euler::Euler( Settings* settings ) : Problem( settings ), _problemType( I_SOD ) {
    _nStates = 3;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( true );
    _sigma = _settings->GetSigma();
    try {
        auto file     = cpptoml::parse_file( _settings->GetInputFile() );
        auto general  = file->get_table( "general" );
        auto ICString = general->get_as<std::string>( "testCase" );
        if( ICString ) {
            if( ICString->compare( "sod" ) == 0 ) {
                _problemType = ICEuler1DType::I_SOD;
            }
            else if( ICString->compare( "highDensityShock" ) == 0 ) {
                _problemType = ICEuler1DType::I_HIGH;
            }
            else {
                _log->error( "[euler1d] Unknown testcase defined!" );
            }
        }

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[Euler] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

Euler::~Euler() {}

void Euler::Solve() {}

Vector Euler::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& /*n*/ ) {
    double rhoInv = 1.0 / u[0];
    double vU     = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( vU, 2 ) );
    double aU     = sqrt( _gamma * p * rhoInv );
    rhoInv        = 1.0 / v[0];
    double vV     = v[1] * rhoInv;
    p             = ( _gamma - 1.0 ) * ( v[2] - 0.5 * v[0] * pow( vV, 2 ) );
    double aV     = sqrt( _gamma * p * rhoInv );

    double uUProjected = nUnit[0] * vU;
    double uVProjected = nUnit[0] * vV;

    double lambdaMin = uUProjected - aU;
    double lambdaMax = uVProjected + aV;

    if( lambdaMin >= 0 )
        return F( u ) * nUnit;
    else if( lambdaMax <= 0 )
        return F( v ) * nUnit;
    else {
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( lambdaMax * F( u ) * nUnit - lambdaMin * F( v ) * nUnit + lambdaMax * lambdaMin * ( v - u ) );
    }
}

Matrix Euler::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix Euler::F( const Vector& u ) {
    double rhoInv = 1.0 / u[0];
    double v      = u[1] * rhoInv;
    double p      = ( _gamma - 1.0 ) * ( u[2] - 0.5 * u[0] * pow( v, 2 ) );
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = u[1] * v + p;
    flux( 2, 0 ) = ( u[2] + p ) * v;
    return flux;
}

Matrix Euler::F( const Matrix& /*u*/ ) {
    _log->error( "[euler] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double Euler::ComputeDt( const Tensor& u, double dx, unsigned level ) const {
    double dtMinTotal = 1e10;
    double dtMin;
    double rhoInv, v, p, a, cfl;
    unsigned kEnd = _settings->GetNqPEAtRef( level );

    cfl = _settings->GetCFL();
    for( unsigned n = 0; n < _settings->GetNMultiElements(); ++n ) {
        for( unsigned k = 0; k < kEnd; ++k ) {
            rhoInv = 1.0 / u( 0, n, k );
            v      = u( 1, n, k ) * rhoInv;
            p      = ( _gamma - 1.0 ) * ( u( 2, n, k ) - 0.5 * u( 0, n, k ) * pow( v, 2 ) );
            a      = sqrt( _gamma * p * rhoInv );

            dtMin      = ( cfl * dx ) * std::min( std::fabs( 1.0 / ( v - a ) ), std::fabs( 1.0 / ( v + a ) ) );
            dtMinTotal = std::min( dtMin, dtMinTotal );
        }
    }

    return dtMinTotal;
}

Vector Euler::IC( const Vector& x, const Vector& xi ) {

    double x0    = 0.5;
    double gamma = 1.4;

    double rhoL = 1.0;
    double rhoR = 0.125;
    double pL   = 1.0;
    double pR   = 0.1;
    double uL   = 0.0;
    double uR   = 0.0;

    if( _problemType == I_HIGH ) {
        rhoR = 0.8;
        pR   = 0.125;    // 0.3;
    }

    Vector y( _nStates );
    if( x[0] < x0 + _sigma[0] * xi[0] ) {
        y[0]                  = rhoL;
        y[1]                  = rhoL * uL;
        double kineticEnergyL = 0.5 * rhoL * pow( uL, 2 );
        double innerEnergyL   = pL / ( gamma - 1.0 );
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

Vector Euler::LoadIC( const Vector& /*x*/, const Vector& /*xi*/ ) {
    _log->error( "[Euler: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}

Matrix Euler::ExactSolution( double t, const Matrix& x, const Vector& xi ) const {
    double x0    = 0.5 + _sigma[0] * xi[0];    // initial shock position
    double rho_l = 1.0;
    double P_l   = 1.0;
    double u_l   = 0.0;
    double rho_r = 0.125;
    double P_r   = 0.1;
    double u_r   = 0.0;

    if( _problemType == I_HIGH ) {
        P_r   = 0.125;    // 0.3
        rho_r = 0.8;
    }

    if( xi.size() > 1 ) {
        rho_r += _sigma[1] * xi[1];
    }
    if( xi.size() > 2 ) {
        P_r += _sigma[2] * xi[2];
    }

    double mu = sqrt( ( _gamma - 1 ) / ( _gamma + 1 ) );

    // speed of sound
    double c_l = sqrt( ( _gamma * P_l / rho_l ) );
    // double c_r = sqrt( ( _gamma * P_r / rho_r ) );

    // find_zero( xi->sod_func( xi, rho_l, P_l, u_l, rho_r, P_r, u_r, _gamma ), ( P_l, P_r ), Bisection() );
    double P_post     = Bisection( P_l, P_r, rho_l, P_l, rho_r, P_r );
    double rho_post   = rho_r * ( ( ( P_post / P_r ) + pow( mu, 2 ) ) / ( 1 + mu * mu * ( P_post / P_r ) ) );
    double rho_middle = rho_l * pow( P_post / P_l, 1 / _gamma );

    // double gm1 = _gamma - 1.0;
    double gp1 = _gamma + 1.0;
    // double gmfac1 = 0.5 * gm1 / _gamma;
    double gmfac2 = 0.5 * gp1 / _gamma;

    double z  = ( P_post / P_r - 1.0 );
    double c5 = sqrt( _gamma * P_r / rho_r );

    double fact = sqrt( 1.0 + gmfac2 * z );

    // shock speed
    double w = c5 * fact;

    double u4 = c5 * z / ( _gamma * fact );

    // Key Positions
    double x4 = x0 + w * t;

    double x1 = x0 - c_l * t;
    double x3 = x0 + u4 * t;
    // determining x2
    double c_2 = c_l - ( ( _gamma - 1.0 ) / 2 ) * u4;
    double x2  = x0 + ( u4 - c_2 ) * t;

    // start setting values
    unsigned nCells = _settings->GetNumCells();

    Vector data_rho( nCells );    // density
    Vector data_P( nCells );      // pressure
    Vector data_u( nCells );      // velocity
    Vector data_e( nCells );      // internal energy

    Matrix fullSolution( nCells, _settings->GetNStates() );

    for( unsigned j = 0; j < nCells; ++j ) {
        if( x( j, 0 ) < x1 ) {
            // Solution b4 x1
            data_rho[j] = rho_l;
            data_P[j]   = P_l;
            data_u[j]   = u_l;
        }
        else if( x1 <= x( j, 0 ) && x( j, 0 ) <= x2 ) {
            // Solution b / w x1 and x2
            double c    = mu * mu * ( ( x0 - x( j, 0 ) ) / t ) + ( 1 - mu * mu ) * c_l;
            data_rho[j] = rho_l * pow( c / c_l, 2 / ( _gamma - 1 ) );
            data_P[j]   = P_l * pow( data_rho[j] / rho_l, _gamma );
            data_u[j]   = ( 1 - mu * mu ) * ( ( -( x0 - x( j, 0 ) ) / t ) + c_l );
        }
        else if( x2 <= x( j, 0 ) && x( j, 0 ) <= x3 ) {
            // Solution b / w x2 and x3
            data_rho[j] = rho_middle;
            data_P[j]   = P_post;
            data_u[j]   = u4;
        }
        else if( x3 <= x( j, 0 ) && x( j, 0 ) <= x4 ) {
            // Solution b / w x3 and x4
            data_rho[j] = rho_post;
            data_P[j]   = P_post;
            data_u[j]   = u4;
        }
        else if( x4 < x( j, 0 ) ) {
            // Solution after x4
            data_rho[j] = rho_r;
            data_P[j]   = P_r;
            data_u[j]   = u_r;
        }
        data_e[j]            = data_P[j] / ( ( _gamma - 1 ) * data_rho[j] );
        fullSolution( j, 0 ) = data_rho[j];
        fullSolution( j, 1 ) = data_rho[j] * data_u[j];
        fullSolution( j, 2 ) = data_e[j];
    }
    return fullSolution;
}

double Euler::SodFunction( double P, double rho_l, double P_l, double rho_r, double P_r ) const {
    double mu = sqrt( ( _gamma - 1.0 ) / ( _gamma + 1.0 ) );
    return ( P - P_r ) * sqrt( ( 1.0 - mu * mu ) * 1.0 / ( rho_r * ( P + mu * mu * P_r ) ) ) -
           ( pow( P_l, ( _gamma - 1.0 / ( 2 * _gamma ) ) ) - pow( P, ( _gamma - 1.0 ) / ( 2 * _gamma ) ) ) *
               ( ( ( 1.0 - mu * mu * mu * mu ) * pow( P_l, 1.0 / _gamma ) * sqrt( 1.0 / ( mu * mu * mu * mu * rho_l ) ) ) );
}

double Euler::Bisection( double PA, double PB, double rho_l, double P_l, double rho_r, double P_r ) const {
    if( SodFunction( PA, rho_l, P_l, rho_r, P_r ) * SodFunction( PB, rho_l, P_l, rho_r, P_r ) >= 0 ) {
        printf( "Incorrect a and b" );
        return -1.0;
    }

    double PC = PA;
    double e  = 1e-10;

    while( std::abs( PB - PA ) >= e ) {
        PC = ( PA + PB ) / 2;
        if( SodFunction( PC, rho_l, P_l, rho_r, P_r ) == 0.0 ) {
            return PC;
        }
        else if( SodFunction( PC, rho_l, P_l, rho_r, P_r ) * SodFunction( PA, rho_l, P_l, rho_r, P_r ) < 0 ) {
            PB = PC;
        }
        else {
            PA = PC;
        }
    }
    return PC;
}
