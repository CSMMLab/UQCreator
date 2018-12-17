#include "shallowwater.h"

ShallowWater::ShallowWater( Settings* settings ) : Problem( settings ) {
    _nStates = 2;
    settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        _g           = 9.81;
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ShallowWater] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ShallowWater::~ShallowWater() {}

void ShallowWater::Solve() {}

Vector ShallowWater::G( const Vector& q_l, const Vector& q_r, const Vector& nUnit, const Vector& n ) {
    double ubar = ( ( q_l[1] / sqrt( q_l[0] ) + q_r[1] / sqrt( q_r[0] ) ) / ( sqrt( q_l[0] ) + sqrt( q_r[0] ) ) );
    double cbar = sqrt( 0.5 * _g * ( q_l[0] + q_r[0] ) );
    double u_r  = q_r[1] / q_r[0];
    double c_r  = sqrt( _g * q_r[0] );
    double u_l  = q_l[1] / q_l[0];
    double c_l  = sqrt( _g * q_l[0] );

    // Compute Einfeldt speeds
    double s1        = ubar + cbar;
    double s2        = ubar - cbar;
    double s3        = u_l + c_l;
    double s4        = u_l - c_l;
    double sMin1     = std::fmin( s1, s2 );
    double sMin2     = std::fmin( s3, s4 );
    double lambdaMin = std::fmin( sMin1, sMin2 );
    s3               = u_r + c_r;
    s4               = u_r - c_r;
    double sMax1     = std::fmax( s1, s2 );
    double sMax2     = std::fmax( s3, s4 );
    double lambdaMax = std::fmax( sMax1, sMax2 );
    Vector s( 2 );
    s[0] = lambdaMin;
    s[1] = lambdaMax;

    // lambdaMin = u_l - lambdaMin;
    // lambdaMax = u_r + lambdaMax;
    // std::cout << lambdaMin << " " << lambdaMax << std::endl;
    /*
            Vector q_hat( 2 );
            q_hat[0] = ( q_r[1] - q_l[1] - lambdaMax * q_r[0] + lambdaMin * q_l[0] ) / ( lambdaMin - lambdaMax );
            q_hat[1] = ( pow( q_r[1], 2 ) / q_r[0] + 0.5 * _g * pow( q_r[0], 2 ) - ( pow( q_l[1], 2 ) / q_l[0] + 0.5 * _g * pow( q_l[0], 2 ) ) -
                         lambdaMax * q_r[1] + lambdaMin * q_l[1] ) /
                       ( lambdaMin - lambdaMax );
            unsigned num_waves = 2;
            unsigned num_eqn   = 2;
            Matrix wave( num_eqn, num_waves );
            wave( 0, 0 ) = q_hat[0] - q_l[0];
            wave( 0, 1 ) = q_r[0] - q_hat[0];
            wave( 1, 0 ) = q_hat[1] - q_l[1];
            wave( 1, 1 ) = q_r[1] - q_hat[1];
            Vector wave1 = q_hat - q_l;
            Vector wave2 = q_r - q_hat;

            Matrix amdq( num_eqn, 1, 0.0 );
            Matrix apdq( num_eqn, 1, 0.0 );

            // Compute variations
            Vector s_index( num_eqn, 0.0 );
            for( unsigned m = 0; m < num_eqn; ++m ) {
                for( unsigned mw = 0; mw < num_waves; ++mw ) {
                    s_index[0] = s[mw];
                    amdq( m, 0 ) += fmin( s_index[0], s_index[1] ) * wave( m, mw );
                    apdq( m, 0 ) += fmax( s_index[0], s_index[1] ) * wave( m, mw );
                }
            }

            return ( apdq + amdq ) * n;*/

    // q [m, LL:UL] -= dtdx [LL:UL] * apdq [m, LL - 1:UL - 1];
    // q [m, LL - 1:UL - 1] -= dtdx [LL - 1:UL - 1] * amdq [m, LL - 1:UL - 1]

    // return F( q_r + lambdaMin * wave1 ) * n;

    double dtdx = _settings->GetCFL() / 12.0;

    return 0.5 * ( F( q_l ) * n + F( q_r ) * n ) - ( 0.5 / dtdx ) * ( q_r - q_l );

    // pow( u[1], 2 ) / u[0] + 0.5 * _g* pow( u[0], 2 );

    // return F( q_hat ) * n;

    lambdaMin = u_l - c_l;
    lambdaMax = u_r + c_r;

    if( lambdaMin >= 0 ) {
        return F( q_l ) * n;
    }
    else if( lambdaMax <= 0 ) {
        return F( q_r ) * n;
    }
    else {
        // return F( ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * q_r - lambdaMin * q_l ) + F( q_l ) * nUnit - F( q_r ) * nUnit ) ) * n;
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( q_l ) * n - lambdaMin * F( q_r ) * n + lambdaMax * lambdaMin * ( q_r - q_l ) * norm( n ) );
    }
}
/*
Vector ShallowWater::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    double uU = u[1] / u[0];
    double uV = v[1] / v[0];
    double cU = sqrt( _g * u[0] );
    double cV = sqrt( _g * v[0] );
    std::cout << cV << " " << cU << std::endl;

    double lambdaMin = std::fmin( uU - cU, uV - cV );
    lambdaMin        = std::fmin( lambdaMin, 0.0 ) - 1.0;
    double lambdaMax = std::fmax( uU + cU, uV + cV );
    lambdaMax        = std::fmax( lambdaMax, 0.0 ) + 0.0;

    // std::cout << lambdaMin << " " << lambdaMax << std::endl;

    // double lambdaMin = uU - cU;
    // double lambdaMax = uV + cV;
    // return F( ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * v - lambdaMin * u ) + F( u ) * nUnit - F( v ) * nUnit ) ) * n;

    if( lambdaMin >= 0 ) {
        return F( u ) * n;
    }
    else if( lambdaMax <= 0 ) {
        return F( v ) * n;
    }
    else if( lambdaMin <= 0 && lambdaMax >= 0 ) {
        // return F( ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * v - lambdaMin * u ) + F( u ) * nUnit - F( v ) * nUnit ) ) * n;
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }
    else {
        exit( EXIT_FAILURE );
    }

    double dtdx = _settings->GetCFL() / 12.0;

    return 0.5 * ( F( u ) * n + F( v ) * n ) - ( 0.5 / dtdx ) * ( v - u );

    double aMinus = std::fmin( uU - cU, uV - cV );
    aMinus        = std::fmin( aMinus, 0.0 );
    double aPlus  = std::fmax( uU + cU, uV + cV );
    aPlus         = std::fmax( aPlus, 0.0 );
    return ( 1.0 / ( aPlus - aMinus ) ) * ( aPlus * F( u ) * n - aMinus * F( v ) * n + aPlus * aMinus * ( v - u ) * norm( n ) );

    // U1(j) = (un1(j - 1) + un1(j + 1))/2 - 0.5 * c * (F1(j+1) - F1(j-1));
    // U2(j) = (un2(j - 1) + un2(j + 1))/2 - 0.5 * c *(F2(j+1) - F2(j-1));
}*/
/*
Vector ShallowWater::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    double hL = u[0];
    double hR = v[0];
    double uL = u[1] / u[0];
    double uR = v[1] / v[0];
    double aL = sqrt( _g * u[0] );
    double aR = sqrt( _g * v[0] );

    double hStar = 0.5 * ( hL + hR ) - 0.25 * ( uR - uL ) * ( hL + hR ) / ( aL + aR );
    double uStar = 0.5 * ( uL + uR ) - ( hR - hL ) * ( aL + aR ) / ( hL + hR );

    double qL, qR;
    if( hStar > hL ) {
        qL = sqrt( 0.5 * ( hStar + hL ) * hStar / pow( hL, 2 ) );
    }
    else
        qL = 1;
    if( hStar > hR ) {
        qR = sqrt( 0.5 * ( hStar + hR ) * hStar / pow( hR, 2 ) );
    }
    else
        qR = 1;
    double sL = uL - aL * qL;
    double sR = uR + aR * qR;

    double lambdaMax = sR;
    double lambdaMin = sL;
    if( lambdaMin >= 0 ) {
        // std::cout << "min" << std::endl;
        return F( u ) * n;
    }
    else if( lambdaMax <= 0 ) {
        // std::cout << "max" << std::endl;
        return F( v ) * n;
    }
    else {
        // return F( ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * v - lambdaMin * u ) + F( u ) * nUnit - F( v ) * nUnit ) ) * n;
        return ( 1.0 / ( lambdaMax - lambdaMin ) ) *
               ( lambdaMax * F( u ) * n - lambdaMin * F( v ) * n + lambdaMax * lambdaMin * ( v - u ) * norm( n ) );
    }

    // Vector ULStar( 2 );
    // double factor = (sL-uL)/(sL-sStar)
    // ULStar[0] =
}*/

Matrix ShallowWater::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix ShallowWater::F( const Vector& u ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = u[1];
    flux( 1, 0 ) = pow( u[1], 2 ) / u[0] + 0.5 * _g * pow( u[0], 2 );
    return flux;
}

Matrix ShallowWater::F( const Matrix& u ) {
    _log->error( "[ShallowWater] Flux not implemented" );
    exit( EXIT_FAILURE );
}

double ShallowWater::ComputeDt( Vector& u, double dx ) const {
    _log->error( "[ShallowWater] ComputeDt not implemented" );
    return 0.0;
}
