#include "pnequations.h"

PNEquations::PNEquations( Settings* settings ) : Problem( settings ), _N( 2 ) {
    _nStates = 4;
    _settings->SetNStates( _nStates );
    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        SetupSystemMatrices();
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[pnequations] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

PNEquations::~PNEquations() {}

double PNEquations::AParam( int l, int k ) const { return std::sqrt( ( l - k + 1 ) * ( l + k + 1 ) / ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ); }

double PNEquations::BParam( int l, int k ) const { return std::sqrt( ( l - k ) * ( l + k ) / ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ); }

double PNEquations::CParam( int l, int k ) const { return std::sqrt( ( l + k + 1 ) * ( l + k + 2 ) / ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ); }

double PNEquations::DParam( int l, int k ) const { return std::sqrt( ( l - k ) * ( l - k - 1 ) / ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ); }

double PNEquations::EParam( int l, int k ) const { return std::sqrt( ( l - k + 1 ) * ( l - k + 2 ) / ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ); }

double PNEquations::FParam( int l, int k ) const { return std::sqrt( ( l + k ) * ( l + k - 1 ) / ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ); }

double PNEquations::CTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * CParam( l, k );
    else
        return CParam( l, k );
}

double PNEquations::DTilde( int l, int k ) const {
    if( k < 0 ) return 0.0;
    if( k == 0 )
        return std::sqrt( 2 ) * DParam( l, k );
    else
        return DParam( l, k );
}

double PNEquations::ETilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * EParam( l, k );
    else
        return EParam( l, k );
}

double PNEquations::FTilde( int l, int k ) const {
    if( k == 1 )
        return std::sqrt( 2 ) * FParam( l, k );
    else
        return FParam( l, k );
}

int PNEquations::Sgn( int k ) const {
    if( k >= 0 )
        return 1;
    else
        return -1;
}

unsigned PNEquations::GlobalIndex( int l, int k ) const {
    unsigned numIndicesPrevLevel  = unsigned( l * l );    // number of previous indices untill level l-1
    unsigned prevIndicesThisLevel = unsigned( k + l );    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

int PNEquations::kPlus( int k ) const { return k + Sgn( k ); }

int PNEquations::kMinus( int k ) const { return k - Sgn( k ); }

void PNEquations::SetupSystemMatrices() {
    unsigned nTotalEntries = GlobalIndex( _N, _N ) + 1;    // total number of entries for sytem matrix
    _Ax                    = Matrix( nTotalEntries, nTotalEntries );
    _Ay                    = Matrix( nTotalEntries, nTotalEntries );
    // loop over columns of A
    for( int l = 0; l < _N; ++l ) {
        for( int k = -l; k <= l; ++k ) {
            // flux matrix in direction x
            if( k != -1 ) {
                _Ax( GlobalIndex( l, k ), GlobalIndex( l - 1, kMinus( k ) ) ) = 0.5 * CTilde( l - 1, std::abs( k ) - 1 );
                _Ax( GlobalIndex( l, k ), GlobalIndex( l + 1, kMinus( k ) ) ) = -0.5 * DTilde( l + 1, std::abs( k ) - 1 );
            }
            _Ax( GlobalIndex( l, k ), GlobalIndex( l - 1, kPlus( k ) ) ) = -0.5 * ETilde( l - 1, std::abs( k ) + 1 );
            _Ax( GlobalIndex( l, k ), GlobalIndex( l + 1, kPlus( k ) ) ) = 0.5 * FTilde( l + 1, std::abs( k ) + 1 );
            // flux matrix in direction y
            if( k != 1 ) {
                _Ay( GlobalIndex( l, k ), GlobalIndex( l - 1, -kMinus( k ) ) ) = -0.5 * CTilde( l - 1, std::abs( k ) - 1 );
                _Ay( GlobalIndex( l, k ), GlobalIndex( l + 1, -kMinus( k ) ) ) = 0.5 * DTilde( l + 1, std::abs( k ) - 1 );
            }
            _Ay( GlobalIndex( l, k ), GlobalIndex( l - 1, -kPlus( k ) ) ) = -0.5 * ETilde( l - 1, std::abs( k ) + 1 );
            _Ay( GlobalIndex( l, k ), GlobalIndex( l + 1, -kPlus( k ) ) ) = 0.5 * FTilde( l + 1, std::abs( k ) + 1 );
        }
    }
}

Vector PNEquations::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {

    Vector uMean = 0.5 * ( u + v );

    return F( uMean ) * n - ( v - u ) * norm( n );
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

Matrix PNEquations::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n ) {
    unsigned nStates = static_cast<unsigned>( u.rows() );
    unsigned Nq      = static_cast<unsigned>( u.columns() );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        column( y, k ) = G( column( u, k ), column( v, k ), nUnit, n );
    }
    return y;
}

Matrix PNEquations::F( const Vector& u ) {
    Matrix flux( u.size(), 2 );

    column( flux, 0 ) = _Ax * u;
    column( flux, 1 ) = _Ay * u;

    return flux;
}

Matrix PNEquations::F( const Matrix& u ) {
    _log->error( "[euler2d] Flux not implemented" );
    exit( EXIT_FAILURE );
    return 0.5 * pow( u, 2 );
}

double PNEquations::ComputeDt( const Matrix& u, double dx ) const {

    double cfl = _settings->GetCFL();

    double maxVelocity = maxVelocity;

    return ( cfl / dx ) * 1 / maxVelocity;
}

Vector PNEquations::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates );
    _sigma            = _settings->GetSigma();
    bool pipeTestCase = true;
    if( pipeTestCase ) {    // pipe testcase
        double gamma = 1.4;
        double R     = 287.87;
        double T     = 273.15;
        double p     = 101325.0;
        double Ma    = 0.0;
        double a     = sqrt( gamma * R * T );

        double uMax  = Ma * a;
        double angle = ( 1.25 + _sigma[0] * 0.0 ) * ( 2.0 * M_PI ) / 360.0;
        double uF    = uMax * cos( angle );
        double vF    = uMax * sin( angle );

        double rhoFarfield = p / ( R * T );

        y[0] = rhoFarfield;
        if( xi.size() > 1 ) {
            y[0] += _sigma[1] * xi[1];
        }
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * y[0] * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( p / ( y[0] * ( gamma - 1 ) ) ) * y[0];
        y[3]                  = kineticEnergyL + innerEnergyL;

        if( x[1] < 1.1 + _sigma[0] * xi[0] ) {
            y[0] = 0.5 * rhoFarfield;
            if( xi.size() > 2 ) y[0] += _sigma[2] * xi[2];
            y[1]           = rhoFarfield * uF;
            y[2]           = rhoFarfield * vF;
            kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
            innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
            y[3]           = 0.5 * ( kineticEnergyL + innerEnergyL );
        }
        /*
         if( x[1] < 1.1 + sigma && x[1] < 1.1 - sigma ) {
             y[0] = 0;
         }
         if( x[1] < 1.1 + sigma && x[1] > 1.1 - sigma ) {
             y[0] = 1;
         }
         if( x[1] > 1.1 + sigma && x[1] > 1.1 - sigma ) {
             y[0] = 2;
         }*/
        return y;
    }
    else {
        double gamma = 1.4;
        double R     = 287.87;
        double T     = 273.15;
        double p     = 101325.0;
        double Ma    = 0.8;
        if( xi.size() == 2 ) {
            Ma = Ma + xi[1] * _sigma[1];
        }
        double a = sqrt( gamma * R * T );

        double uMax  = Ma * a;
        double angle = ( 1.25 + _sigma[0] * xi[0] ) * ( 2.0 * M_PI ) / 360.0;
        double uF    = uMax * cos( angle );
        double vF    = uMax * sin( angle );

        double rhoFarfield = p / ( R * T );

        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        double kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        double innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        return y;
    }
}
