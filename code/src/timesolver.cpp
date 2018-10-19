#include "timesolver.h"
#include "expliciteuler.h"

TimeSolver::TimeSolver( Settings* settings, Mesh* mesh ) : _settings( settings ), _mesh( mesh ) {
    _CFL        = _settings->GetCFL();
    _dx         = _mesh->GetArea( 0 );
    _dt         = _dx * _settings->GetCFL() / 12.0;
    _nTimeSteps = static_cast<unsigned>( _settings->GetTEnd() / _dt );
}

TimeSolver::~TimeSolver() {}

TimeSolver* TimeSolver::Create( Settings* settings, Mesh* mesh ) {
    auto file          = cpptoml::parse_file( settings->GetInputFile() );
    auto section       = file->get_table( "problem" );
    std::string method = section->get_as<std::string>( "timestepping" ).value_or( "" );
    if( method.compare( "explicitEuler" ) == 0 || method.compare( "EE" ) == 0 ) {
        return new ExplicitEuler( settings, mesh );
    }
    else {
        std::cerr << "Invalid timesolver type" << std::endl;
        exit( EXIT_FAILURE );
        return nullptr;
    }
}

double TimeSolver::GetTimeStepSize() { return _dt; }

double TimeSolver::GetNTimeSteps() { return _nTimeSteps; }
