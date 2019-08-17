#include "quadraturegrid.h"

#include "tensorizedcc.h"
#include "tensorizedquadrature.h"
#include "uniformsparsegrid.h"

QuadratureGrid::QuadratureGrid() {}

QuadratureGrid* QuadratureGrid::Create( Settings* settings ) {
    auto log      = spdlog::get( "event" );
    auto gridType = settings->GetGridType();
    if( gridType == GridType::G_SPARSEGRID ) {
        return new UniformSparseGrid( settings );
    }
    else if( gridType == GridType::G_TENSORIZEDGRID ) {
        return new TensorizedQuadrature( settings );
    }
    else if( gridType == GridType::G_TENSORIZEDCC ) {
        return new TensorizedCC( settings );
    }
    else {
        log->error( "[quadratureGrid]: Invalid grid type" );
        exit( EXIT_FAILURE );
    }
}

QuadratureGrid* QuadratureGrid::Create( Settings* settings, unsigned level, unsigned dim ) {
    auto log      = spdlog::get( "event" );
    auto gridType = settings->GetGridType();
    if( gridType == GridType::G_SPARSEGRID ) {
        return new UniformSparseGrid( level, dim );
    }
    else if( gridType == GridType::G_TENSORIZEDGRID ) {
        std::cerr << "G_TENSORIZEDGRID Not Implemented" << std::endl;
        return new TensorizedQuadrature( settings, level, dim );
    }
    else if( gridType == GridType::G_TENSORIZEDCC ) {
        std::cerr << "G_TENSORIZEDCC Not Implemented" << std::endl;
        return new TensorizedCC( settings, level );
    }
    else {
        log->error( "[quadratureGrid]: Invalid grid type" );
        exit( EXIT_FAILURE );
    }
}
