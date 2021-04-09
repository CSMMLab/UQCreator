#include "problems/externalproblem.h"

ExternalProblem::ExternalProblem( Settings* settings ) : Problem( settings ), _problemType( I_NACA ) {
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
                _log->error( "[ExternalProblem] Unknown testcase defined!" );
            }
        }

        auto problem = file->get_table( "problem" );
        _gamma       = problem->get_as<double>( "gamma" ).value_or( 1.4 );
        _settings->SetGamma( _gamma );
    } catch( const cpptoml::parse_exception& e ) {
        _log->error( "[ExternalProblem] Failed to parse {0}: {1}", _settings->GetInputFile(), e.what() );
        exit( EXIT_FAILURE );
    }
}

ExternalProblem::~ExternalProblem() {}

void ExternalProblem::Solve() {}

Vector ExternalProblem::G( const Vector& u, const Vector& v, const Vector& nUnit, const Vector& n ) {
}

Matrix ExternalProblem::G( const Matrix& u, const Matrix& v, const Vector& nUnit, const Vector& n, unsigned level ) {
}

Matrix ExternalProblem::F( const Vector& u ) const {
}

Matrix ExternalProblem::F( const Matrix& u ) {
}

Matrix ExternalProblem::BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const {
}

double ExternalProblem::ComputeDt( const Tensor& u, double dx, unsigned level ) const {
}

Vector ExternalProblem::IC( const Vector& x, const Vector& xi ) {
}

Vector ExternalProblem::LoadIC( const Vector& x, const Vector& xi ) {
}
