#include "polynomial.h"
#include "hermite.h"
#include "legendre.h"

Polynomial::Polynomial( unsigned degree ) : _degree( degree ) {}

void Polynomial::Sort() {
    std::vector<std::size_t> p( _nodes.size() );
    std::iota( p.begin(), p.end(), 0 );
    std::sort( p.begin(), p.end(), [&]( std::size_t i, std::size_t j ) { return _nodes[i] < _nodes[j]; } );
    Vector sorted_nodes( p.size() ), sorted_weights( p.size() );
    std::transform( p.begin(), p.end(), sorted_nodes.begin(), [&]( std::size_t i ) { return _nodes[i]; } );
    std::transform( p.begin(), p.end(), sorted_weights.begin(), [&]( std::size_t i ) { return _weights[i]; } );
    _nodes   = sorted_nodes;
    _weights = sorted_weights;
}

Polynomial* Polynomial::Create( Settings* settings, unsigned order ) {
    auto log = spdlog::get( "event" );
    if( settings->GetDistributionType() == DistributionType::D_LEGENDRE ) {
        return new Legendre( order );
    }
    else if( settings->GetDistributionType() == DistributionType::D_HERMITE ) {
        return new Hermite( settings->GetNQuadPoints() );
    }
    else {
        log->error( "[distribution] Invalid distribution type" );
        exit( EXIT_FAILURE );
    }
}
