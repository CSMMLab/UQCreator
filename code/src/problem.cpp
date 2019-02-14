#include <spdlog/spdlog.h>

#include "burgers.h"
#include "euler.h"
#include "euler2d.h"
#include "problem.h"
#include "shallowwater.h"
#include "shallowwater2d.h"

Problem::Problem( Settings* settings ) : _settings( settings ) {
    _log = spdlog::get( "event" );
    _settings->SetExactSolution( false );
}

Problem::~Problem() {}

Problem* Problem::Create( Settings* settings ) {
    auto log = spdlog::get( "event" );
    if( settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
        return new Burgers( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_EULER_1D ) {
        return new Euler( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_EULER_2D ) {
        return new Euler2D( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_SHALLOWWATER_1D ) {
        return new ShallowWater( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_SHALLOWWATER_2D ) {
        return new ShallowWater2D( settings );
    }
    else {
        log->error( "[Problem] Invalid problem type!" );
        exit( EXIT_FAILURE );
    }
}

double Problem::ComputeDt( const Matrix& u, double dx ) const { return _settings->GetDT(); }

Vector Problem::ExactSolution( double t, const Vector& x, const Vector& xi ) const {
    _log->error( "[Problem]: No exact solution specified" );
    exit( EXIT_FAILURE );
}
