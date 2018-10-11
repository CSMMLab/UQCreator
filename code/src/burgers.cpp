#include "burgers.h"

Burgers::Burgers( std::string inputFile ) : Problem( inputFile ) {
    _x          = _mesh->GetGrid();
    _dx         = _mesh->GetSpacing().at( 0 );
    _dt         = _dx * _CFL / 12.0;
    _nCells     = _mesh->GetNumCells();
    _nStates    = 1;
    _nTimeSteps = static_cast<unsigned>( _tEnd / _dt );
    _tEnd       = _nTimeSteps * _dt;
    _u          = Vector( _nCells + 4, 0.0 );
}

double Burgers::H( double u, double v, double w ) { return v - ( _dt / _dx ) * ( G( v, w ) - G( u, v ) ); }

double Burgers::G( double u, double v ) { return F( u ); }

Matrix Burgers::G( const Matrix& u, const Matrix& v ) { return F( u ); }

double Burgers::F( double u ) { return 0.5 * u * u; }

Matrix Burgers::F( const Matrix& u ) { return 0.5 * blaze::pow( u, 2 ); }

void Burgers::Solve() {
    double* uNew = new double[_nCells + 4];

    // setup IC
    for( unsigned j = 0; j < _nCells + 4; ++j ) {
        _u[j] = IC( _x[j], 12.0, 3.0 );
    }

    for( unsigned n = 0; n < _nTimeSteps; ++n ) {
        for( unsigned j = 2; j < _nCells + 2; ++j ) {
            //_dt = _timeDiscretization->getDt();
            uNew[j] = H( _u[j - 1], _u[j], _u[j + 1] );
        }

        for( unsigned j = 2; j < _nCells + 2; ++j ) {
            _u[j] = uNew[j];
        }
    }
}

double Burgers::IC( double x, double uL, double uR ) {
    double a = 0.5;
    double b = 1.5;
    if( x < a ) {
        return uL;
    }
    else if( x > a && x < b ) {
        return uR + ( uL - uR ) * ( b - x ) / ( b - a );
    }
    else {
        return uR;
    }
}

void Burgers::Print() {
    std::ofstream out( "outFile" );
    for( unsigned j = 2; j < _nCells + 2; ++j ) {
        out << _x[j] << " " << _u[j] << std::endl;
    }
}

void Burgers::WriteToFile( std::string filename, int filetype ) const {}

double Burgers::ExactSolution( double t, double x, double xi ) {
    double uL    = 12.0;
    double uR    = 3.0;
    double sigma = 0.2;
    double x0    = 0.5;
    double x1    = 1.5;
    double x0Bar;
    double x1Bar;
    bool shock = true;
    double y   = 0;
    if( shock ) {
        x0Bar = x0 + sigma * xi + uL * ( 1.0 / 9.0 ) + 0.5 * ( uL + uR ) * ( t - ( 1.0 / 9.0 ) );
        x1Bar = x1 + sigma * xi + uR * ( 1.0 / 9.0 ) + 0.5 * ( uL + uR ) * ( t - ( 1.0 / 9.0 ) );
        if( x <= x0Bar ) {
            return uL;
        }
        else if( x > x1Bar ) {
            return uR;
        }
    }
    else {
        x0Bar = x0 + sigma * xi + uL * t;
        x1Bar = x1 + sigma * xi + uR * t;
        if( x < x0Bar ) {
            y = uL;
        }
        else if( x > x1Bar ) {
            y = uR;
        }
        else {
            y = uL + ( uR - uL ) * ( x0Bar - x ) / ( x0Bar - x1Bar );
        }
    }
    return y;
}
