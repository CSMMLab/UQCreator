#include "pnequations.h"

PNEquations::PNEquations( Settings* settings ) : Problem( settings ), _N( 7 ) {
    _nStates = unsigned( GlobalIndex( _N, _N ) + 1 );
    _settings->SetNStates( _nStates );
    _settings->SetSource( true );
    _sigmaA = 0.0;    // absorption coefficient
    _sigmaS = 1.0;    // scattering coefficient
    _sigmaT = _sigmaA + _sigmaS;
    try {
        auto file    = cpptoml::parse_file( _settings->GetInputFile() );
        auto problem = file->get_table( "problem" );
        SetupSystemMatrices();
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[pnequations] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

PNEquations::PNEquations( Settings* settings, bool noSystemMatrix ) : Problem( settings ) {}

PNEquations::~PNEquations() {}

double PNEquations::AParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l + k + 1 ) ) / double( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) );
}

double PNEquations::BParam( int l, int k ) const {
    return std::sqrt( double( ( l - k ) * ( l + k ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
}

double PNEquations::CParam( int l, int k ) const {
    return std::sqrt( double( ( l + k + 1 ) * ( l + k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double PNEquations::DParam( int l, int k ) const {
    return std::sqrt( double( ( l - k ) * ( l - k - 1 ) ) / double( ( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) ) );
}

double PNEquations::EParam( int l, int k ) const {
    return std::sqrt( double( ( l - k + 1 ) * ( l - k + 2 ) ) / double( ( ( 2 * l + 3 ) * ( 2 * l + 1 ) ) ) );
}

double PNEquations::FParam( int l, int k ) const {
    return std::sqrt( double( ( l + k ) * ( l + k - 1 ) ) / double( ( 2 * l + 1 ) * ( 2 * l - 1 ) ) );
}

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
    else {
        return FParam( l, k );
    }
}

int PNEquations::Sgn( int k ) const {
    if( k >= 0 )
        return 1;
    else
        return -1;
}

int PNEquations::GlobalIndex( int l, int k ) const {
    int numIndicesPrevLevel  = l * l;    // number of previous indices untill level l-1
    int prevIndicesThisLevel = k + l;    // number of previous indices in current level
    return numIndicesPrevLevel + prevIndicesThisLevel;
}

int PNEquations::kPlus( int k ) const { return k + Sgn( k ); }

int PNEquations::kMinus( int k ) const { return k - Sgn( k ); }

void PNEquations::SetupSystemMatrices() {
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
                if( j >= 0 && j < int( nTotalEntries ) ) _Ax( i, unsigned( j ) ) = 0.5 * CTilde( l - 1, std::abs( k ) - 1 );

                j = GlobalIndex( l + 1, kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) ) _Ax( i, unsigned( j ) ) = -0.5 * DTilde( l + 1, std::abs( k ) - 1 );
            }

            j = GlobalIndex( l - 1, kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) ) _Ax( i, unsigned( j ) ) = -0.5 * ETilde( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) ) _Ax( i, unsigned( j ) ) = 0.5 * FTilde( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction y
            if( k != 1 ) {
                j = GlobalIndex( l - 1, -kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) ) _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * CTilde( l - 1, std::abs( k ) - 1 );

                j = GlobalIndex( l + 1, -kMinus( k ) );
                if( j >= 0 && j < int( nTotalEntries ) ) _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * DTilde( l + 1, std::abs( k ) - 1 );
            }

            j = GlobalIndex( l - 1, -kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) ) _Ay( i, unsigned( j ) ) = -0.5 * Sgn( k ) * ETilde( l - 1, std::abs( k ) + 1 );

            j = GlobalIndex( l + 1, -kPlus( k ) );
            if( j >= 0 && j < int( nTotalEntries ) ) _Ay( i, unsigned( j ) ) = 0.5 * Sgn( k ) * FTilde( l + 1, std::abs( k ) + 1 );

            //
            // flux matrix in direction z
            j = GlobalIndex( l - 1, k );
            if( j >= 0 && j < int( nTotalEntries ) ) _Az( i, unsigned( j ) ) = AParam( l - 1, k );

            j = GlobalIndex( l + 1, k );
            if( j >= 0 && j < int( nTotalEntries ) ) _Az( i, unsigned( j ) ) = BParam( l + 1, k );
        }
    }
    std::cout << "System Matrix Set UP!" << std::endl;
    std::cout << "A_x =" << _Ax << std::endl;
    std::cout << "A_y =" << _Ay << std::endl;
    std::cout << "A_z =" << _Az << std::endl;
}

Vector PNEquations::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    return F( 0.5 * ( u + v ) ) * n - 0.5 * ( v - u ) * norm( n );
    // F( 0.5 * ( u + v ) ) * n - 0.5 * ( ( v - u ) * norm( n ) * nUnit[0] + ( v - u ) * norm( n ) * nUnit[1] );
    // return F( 0.5 * ( u + v ) ) * n - 0.5 * ( v - u ) * norm( n ) / _settings->GetDT(); // LF does not work since norm(n) != Area(j)
}

Matrix PNEquations::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
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
    // column( flux, 1 ) = _Az * u;

    return flux;
}

Matrix PNEquations::F( const Matrix& u ) {
    _log->error( "[euler2d] Flux not implemented" );
    exit( EXIT_FAILURE );
    return 0.5 * pow( u, 2 );
}

Matrix PNEquations::Source( const Matrix& uQ, const Vector& x, double t, unsigned level ) const {
    unsigned nStates = static_cast<unsigned>( uQ.rows() );
    unsigned Nq      = static_cast<unsigned>( uQ.columns() );
    Vector g( nStates, 0.0 );
    g[0] = 1.0;    // 2 * M_PI;
    Matrix y( nStates, Nq, 0.0 );
    for( unsigned k = 0; k < Nq; ++k ) {
        for( unsigned s = 0; s < nStates; ++s ) {
            y( s, k ) = -_sigmaA * uQ( s, k ) - _sigmaS * ( 1.0 - g[s] ) * uQ( s, k );
            // y( s, k ) = _sigmaT * uQ( s, k ) - _sigmaS * g[s] * uQ( s, k );
        }
    }
    return y;
}

double PNEquations::ComputeDt( const Matrix& u, double dx, unsigned level ) const {

    double cfl = _settings->GetCFL();

    double maxVelocity = 1.0;

    return ( cfl / dx ) / maxVelocity;
}

Vector PNEquations::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double x0    = 0.0;
    double y0    = 0.0;
    double s2    = 3.2 * std::pow( 0.01, 2 );    // std::pow( 0.03, 2 );
    double floor = 0.0;
    _sigma       = _settings->GetSigma();

    y[0] = std::fmax( floor, 1.0 / ( 4.0 * M_PI * s2 ) * exp( -( ( x[0] - x0 ) * ( x[0] - x0 ) + ( x[1] - y0 ) * ( x[1] - y0 ) ) / 4.0 / s2 ) );

    return y;
}

Vector PNEquations::LoadIC( const Vector& x, const Vector& xi ) { return Vector( 1 ); }
