#include "problem.h"
#include "burgers.h"
#include "euler.h"
#include "euler2d.h"

Problem::Problem( Settings* settings ) : _settings( settings ) {}

Problem::~Problem() {}

Problem* Problem::Create( Settings* settings ) {
    if( settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
        return new Burgers( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_EULER_1D ) {
        return new Euler( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_EULER_2D ) {
        return new Euler2D( settings );
    }
    else {
        std::cerr << "[Problem] Invalid problem type!" << std::endl;
        exit( EXIT_FAILURE );
    }
}
