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
    _settings->SetSource( true );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

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
                if( j >= 0 && j < int( nTotalEntries ) )
                    _Ax( i, unsigned( j ) ) = -0.5 * Delta( l + 1, std::abs( k ) - 1 ) * DTilde( l + 1, std::abs( k ) - 1 );
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
            for( unsigned j = 0; j < nTotalEntries; ++j ) {
                _Ax( i, j ) = _c * _Ax( i, j ) * Delta( l, k );
                _Ay( i, j ) = _c * _Ay( i, j ) * Delta( l, k );
                _Az( i, j ) = _c * _Az( i, j ) * Delta( l, k );
            }
        }
    }
}

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
