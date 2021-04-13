#include "polynomial.h"
#include "hermite.h"
#include "legendre.h"

Polynomial::Polynomial( unsigned degree ) : _degree( degree ) {}

void Polynomial::Sort() {
    std::vector<unsigned> p( _nodes.size() );
    std::iota( p.begin(), p.end(), 0 );
    std::sort( p.begin(), p.end(), [&]( unsigned i, unsigned j ) { return _nodes[i] < _nodes[j]; } );
    Vector sorted_nodes( static_cast<unsigned>( p.size() ) ), sorted_weights( static_cast<unsigned>( p.size() ) );
    std::transform( p.begin(), p.end(), sorted_nodes.begin(), [&]( unsigned i ) { return _nodes[i]; } );
    std::transform( p.begin(), p.end(), sorted_weights.begin(), [&]( unsigned i ) { return _weights[i]; } );
    _nodes   = sorted_nodes;
    _weights = sorted_weights;
}

Polynomial* Polynomial::Create( Settings* /*settings*/, unsigned order, DistributionType distributionType ) {
    auto log = spdlog::get( "event" );
    if( distributionType == DistributionType::D_LEGENDRE ) {
        return new Legendre( order );
    }
    else if( distributionType == DistributionType::D_HERMITE ) {
        return new Hermite( order );
    }
    else {
        log->error( "[distribution] Invalid distribution type" );
        exit( EXIT_FAILURE );
    }
}
