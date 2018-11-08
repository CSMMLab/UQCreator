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
    auto log = spdlog::get( "event" );
    if( settings->GetTimesteppingType() == TimesteppingType::T_EXPLICITEULER ) {
        return new ExplicitEuler( settings, mesh );
    }
    else {
        log->error( "[timesolver] Invalid timesolver type" );
        exit( EXIT_FAILURE );
    }
}

double TimeSolver::GetTimeStepSize() { return _dt; }

double TimeSolver::GetNTimeSteps() { return _nTimeSteps; }
