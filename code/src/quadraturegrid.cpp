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
    else {
        log->error( "[quadratureGrid]: Invalid grid type" );
        exit( EXIT_FAILURE );
    }
}

QuadratureGrid* QuadratureGrid::Create( Settings* settings, unsigned level ) {
    auto log      = spdlog::get( "event" );
    auto gridType = settings->GetGridType();
    if( gridType == GridType::G_SPARSEGRID ) {
        return new UniformSparseGrid( settings, level );
    }
    else if( gridType == GridType::G_TENSORIZEDGRID ) {
        return new TensorizedQuadrature( settings, level );
    }
    else if( gridType == GridType::G_TENSORIZEDCC ) {
        return new TensorizedCC( settings, level );
    }
    else {
        log->error( "[quadratureGrid]: Invalid grid type" );
        exit( EXIT_FAILURE );
    }
}
