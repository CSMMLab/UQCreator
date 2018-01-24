#include "timesolver.h"
#include "thetamethod.h"

TimeSolver::TimeSolver( Problem* problem ) : _problem( problem ) {
    _CFL        = _problem->GetCFL();
    _dx         = _problem->GetMesh()->GetSpacing()[0];
    _dt         = _dx * _problem->GetCFL() / 12.0;
    _nTimeSteps = _problem->GetTEnd() / _dt;
}

TimeSolver::~TimeSolver() {}

TimeSolver* TimeSolver::Create( Problem* problem ) {
    auto file          = cpptoml::parse_file( problem->GetInputFile() );
    auto section       = file->get_table( "problem" );
    std::string method = section->get_as<std::string>( "timestepping" ).value_or( "" );
    if( method.compare( "explicitEuler" ) == 0 ) {
        return new ThetaMethod( problem, 0.0 );
    }
    else {
        std::cerr << "Invalid timesolver type" << std::endl;
        exit( EXIT_FAILURE );
        return NULL;
    }
}

double TimeSolver::GetTimeStepSize() { return _dt; }
double TimeSolver::GetNTimeSteps() { return _nTimeSteps; }
