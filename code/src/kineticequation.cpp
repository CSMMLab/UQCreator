#include "kineticequation.h"
#include "quadraturegrid.h"

KineticEquation::KineticEquation( Settings* settings ) : Problem( settings ) {
    _nStates = 1;
    settings->SetNStates( _nStates );
    _settings->SetExactSolution( false );

    // get quadrature grid
    _grid   = QuadratureGrid::Create( _settings, _settings->GetNQuadPoints() );
    _xiQuad = _grid->GetNodes();

    try {
        auto file = cpptoml::parse_file( _settings->GetInputFile() );

        auto problem = file->get_table( "problem" );
        //_gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        //_settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[KineticEquation] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

KineticEquation::~KineticEquation() { delete _grid; }

Vector KineticEquation::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
    unused( u );
    unused( v );
    unused( nUnit );
    unused( n );

    _log->error( "[KineticEquation: G scalar not implemented]" );
    exit( EXIT_FAILURE );
}

Matrix KineticEquation::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
    unused( n );

    unsigned nStates = u.rows();
    unsigned Nq      = _settings->GetNqPEAtRef( level );
    Matrix y( nStates, Nq );
    for( unsigned k = 0; k < Nq; ++k ) {
        if( _xiQuad[k][0] > 0.0 ) {
            column( y, k ) = F( column( u, k ), _xiQuad[k][0] ) * nUnit;
        }
        else {
            column( y, k ) = F( column( v, k ), _xiQuad[k][0] ) * nUnit;
        }
    }
    return y;
}

Matrix KineticEquation::F( const Vector& u, double mu ) {
    Matrix flux( u.size(), 1 );
    flux( 0, 0 ) = mu * u[0];
    return flux;
}

double KineticEquation::ComputeDt( const Matrix& u, double dx, unsigned level ) const {
    unused( u );
    unused( level );

    double cfl = _settings->GetCFL();
    return cfl * dx;
}

Vector KineticEquation::IC( const Vector& x, const Vector& xi ) {
    Vector y( _nStates, 0.0 );
    double x0    = 0.0;
    double floor = 1e-4;
    _sigma       = _settings->GetSigma();
    double sigmaXi;
    if( xi.size() > 1 )
        sigmaXi = _sigma[1] * xi[1];
    else
        sigmaXi = 0.0;

    y[0] = std::fmax(
        floor, pow( 50.0, 2 ) / ( 8.0 * M_PI * pow( sigmaXi + 2.0, 2 ) ) * exp( -0.5 * pow( 50.0 * ( x[0] - x0 ), 2 ) / pow( sigmaXi + 2.0, 2 ) ) );
    return y;
}

Vector KineticEquation::LoadIC( const Vector& x, const Vector& xi ) {
    unused( x );
    unused( xi );

    _log->error( "[KineticEquation: LoadIC not implemented]" );
    exit( EXIT_FAILURE );
}

Matrix KineticEquation::ExactSolution( double t, const Matrix& x, const Vector& xi ) const {
    unused( t );
    unused( x );
    unused( xi );

    _log->error( "[KineticEquation: ExactSolution not implemented]" );
    exit( EXIT_FAILURE );
}
