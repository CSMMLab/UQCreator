#include "closure.h"
#include "mathtools.h"
#include "stochasticgalerkin.h"
#include "tensorizedquadrature.h"
#include "uniformsparsegrid.h"

Closure::Closure( Settings* settings, QuadratureGrid* quadGrid )
    : _settings( settings ), _nMoments( _settings->GetNMoments() ), _nQuadPoints( _settings->GetNQuadPoints() ), _nStates( _settings->GetNStates() ) {
    _log = spdlog::get( "event" );
    // initialize classes: Basis Functions and Quadrature rules are defined for Legendre and Hermite
    _quad.resize( 2 );
    _quad[0] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_LEGENDRE );
    _quad[1] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_HERMITE );

    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();

    // setup map from k (0,...,_nQTotal-1) to individual indices
    std::vector<std::vector<unsigned>> indices;

    // setup map from i (0,...,nTotal-1) to individual indices
    std::vector<unsigned> indexTest;
    indexTest.resize( _numDimXi );
    unsigned totalDegree;    // compute total degree of basis function i
    for( unsigned i = 0; i < std::pow( _settings->GetNMoments(), _numDimXi ); ++i ) {
        totalDegree = 0;
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indexTest[l] = unsigned( ( i - i % unsigned( std::pow( _nMoments, l ) ) ) / unsigned( std::pow( _nMoments, l ) ) ) % _nMoments;
            totalDegree += indexTest[l];
        }
        // if total degree is sufficiently small or max degree is used, indices are stored
        if( totalDegree < _nMoments || _settings->UsesMaxDegree() ) indices.push_back( indexTest );
    }
    _nTotal = unsigned( indices.size() );

    // TensorizedQuadrature* quadGrid = new TensorizedQuadrature( _settings );
    // std::cout << "quadGrid done" << std::endl;
    auto xiGrid = quadGrid->GetNodes();
    auto wGrid  = quadGrid->GetWeights();

    // set total number of quadrature points
    _nQTotal = quadGrid->GetNodeCount();

    // compute basis functions evaluated at the quadrature points
    _phiTildeWf = Matrix( _nQTotal, _nTotal, 1.0 );

    unsigned n = 0;
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned l = 0; l < _numDimXi; ++l ) {
                if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                _phiTildeWf( k, i ) *=
                    _quad[n]->Evaluate( indices[i][l], xiGrid[k][l] ) / _quad[n]->L2Norm( indices[i][l] ) * _quad[n]->fXi( xiGrid[k][l] );
            }

            _phiTildeWf( k, i ) *= wGrid[k];
        }
    }

    /*
        // Test if polynomials are orthonormal
        Matrix testQuad( _nTotal, _nTotal, 0.0 );
        for( unsigned i = 0; i < _nTotal; ++i ) {
            for( unsigned j = 0; j < _nTotal; ++j ) {
                for( unsigned k = 0; k < _nQTotal; ++k ) {
                    testQuad( i, j ) += _phiTilde( k, i ) * _phiTildeWf( k, j );
                }
            }
        }
        std::cout << "test I " << testQuad << std::endl;*/
}

Closure::~Closure() {
    for( unsigned l = 0; l < _quad.size(); ++l ) {
        delete _quad[l];
    }
}

Closure* Closure::Create( Settings* settings, QuadratureGrid* quadGrid ) {
    auto log         = spdlog::get( "event" );
    auto closureType = settings->GetClosureType();
    if( closureType == ClosureType::C_STOCHASTICGALERKIN ) {
        return new StochasticGalerkin( settings, quadGrid );
    }
    else {
        log->error( "[closure]: Invalid closure type" );
        exit( EXIT_FAILURE );
    }
}
