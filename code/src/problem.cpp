#include <spdlog/spdlog.h>

#include "advection2d.h"
#include "burgers.h"
#include "euler.h"
#include "euler2d.h"
#include "kineticequation.h"
#include "m1equations1d.h"
#include "navierstokes.h"
#include "pnequations.h"
#include "pnequations1d.h"
#include "problem.h"
#include "radiationhydrodynamics.h"
#include "radiationhydrodynamics1d.h"
#include "shallowwater.h"
#include "shallowwater2d.h"
#include "thermalpn.h"
#include "thermalpn2d.h"
#include "thermalradiative.h"
#include "thermalradiativegeneral.h"

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
    else if( settings->GetProblemType() == ProblemType::P_ADVECTION_2D ) {
        return new Advection2D( settings );
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
    else if( settings->GetProblemType() == ProblemType::P_PNEQUATIONS_2D ) {
        return new PNEquations( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_PNEQUATIONS_1D ) {
        return new PNEquations1D( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_M1EQUATIONS_1D ) {
        return new M1Equations1D( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_RADIATIONHYDRO_2D ) {
        return new RadiationHydrodynamics( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_RADIATIONHYDRO_1D ) {
        return new RadiationHydrodynamics1D( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_THERMALRAD_1D ) {
        return new ThermalRadiativeGeneral( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_THERMALPN_1D ) {
        return new ThermalPN( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_THERMALPN_2D ) {
        return new ThermalPN2D( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_KINETIC_1D ) {
        return new KineticEquation( settings );
    }
    else if( settings->GetProblemType() == ProblemType::P_NAVIERSTOKES_1D ) {
        return new NavierStokes( settings );
    }
    else {
        log->error( "[Problem] Invalid problem type!" );
        exit( EXIT_FAILURE );
    }
}

double Problem::ComputeDt( const Tensor& u, double dx, unsigned level ) const {
    unused( u );
    unused( dx );
    unused( level );

    std::cout << "Using dt computation from problem class" << std::endl;
    return _settings->GetDT();
}

Matrix Problem::ExactSolution( double t, const Matrix& x, const Vector& xi ) const {
    unused( t );
    unused( x );
    unused( xi );

    _log->error( "[Problem]: No exact solution specified" );
    exit( EXIT_FAILURE );
}

Matrix Problem::Source( const Matrix& uQ ) const {
    std::cerr << "[Problem]: Source not defined" << std::endl;
    exit( EXIT_FAILURE );
    return uQ;
}

Tensor Problem::Source( const Tensor& uQ, const Vector& x, double t, unsigned level ) const {
    unused( x );
    unused( t );
    unused( level );

    std::cerr << "[Problem]: Source not defined" << std::endl;
    exit( EXIT_FAILURE );
    return uQ;
}

void Problem::SourceImplicit( Matrix& uQNew, const Matrix& uQTilde, const Matrix& uQ, const Vector& x, double t, unsigned level ) const {
    unused( uQNew );
    unused( uQTilde );
    unused( uQ );
    unused( x );
    unused( t );
    unused( level );

    std::cerr << "[Problem]: Source not defined" << std::endl;
    exit( EXIT_FAILURE );
}

Matrix Problem::BoundaryFlux( const Matrix& u, const Vector& nUnit, const Vector& n, unsigned level ) const {
    unused( u );
    unused( nUnit );
    unused( n );
    unused( level );

    std::cerr << "[Problem]: Boundary Flux not defined" << std::endl;
    exit( EXIT_FAILURE );
}
