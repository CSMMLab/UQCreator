#include "tensorizedquadrature.h"

#include "legendre.h"

TensorizedQuadrature::TensorizedQuadrature( Settings* settings ) : _settings( settings ), _nQuadPoints( settings->GetNQuadPoints() ) {
    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();
    _nQTotal  = unsigned( std::pow( _settings->GetNQuadPoints(), _numDimXi ) );
    CreateGrid();
}

TensorizedQuadrature::TensorizedQuadrature( Settings* settings, unsigned nQuadPoints ) : _settings( settings ), _nQuadPoints( nQuadPoints ) {
    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();
    _nQTotal  = unsigned( std::pow( nQuadPoints, _numDimXi ) );
    CreateGrid();
}

TensorizedQuadrature::~TensorizedQuadrature() {
    for( unsigned i = 0; i < _numDimXi; i++ ) delete[] _nodes[i];
    delete[] _nodes;

    delete[] _weights;
}

void TensorizedQuadrature::CreateGrid() {
    // allocate nodes and weights
    _nodes = new double*[_numDimXi];
    for( unsigned i = 0; i < _numDimXi; i++ ) {
        _nodes[i] = new double[_nQTotal];
    }
    for( unsigned i = 0; i < _numDimXi; i++ )
        for( unsigned j = 0; j < _nQTotal; j++ ) _nodes[i][j] = 0.0;

    _weights = new double[_nQTotal];
    for( unsigned i = 0; i < _nQTotal; i++ ) _weights[i] = 1.0;

    _quad.resize( 2 );

    _quad[0] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_LEGENDRE );
    _quad[1] = Polynomial::Create( _settings, _nQuadPoints, DistributionType::D_HERMITE );

    // setup map from k (0,...,_nQTotal-1) to individual indices
    std::vector<std::vector<unsigned>> indicesQ;
    indicesQ.resize( _nQTotal );
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        indicesQ[k].resize( _numDimXi );
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            indicesQ[k][l] = unsigned( ( k - k % unsigned( std::pow( _nQuadPoints, l ) ) ) / unsigned( std::pow( _nQuadPoints, l ) ) ) % _nQuadPoints;
        }
    }

    // store quad points and weights for individual distributions
    std::vector<Vector> xi;
    xi.resize( _quad.size() );
    std::vector<Vector> w;
    w.resize( _quad.size() );
    for( unsigned l = 0; l < _quad.size(); ++l ) {
        xi[l] = _quad[l]->GetNodes();
        w[l]  = _quad[l]->GetWeights();
    }

    // setup weights and nodes
    unsigned n = 0;
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            if( _settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
            if( _settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
            _nodes[l][k] = xi[n][indicesQ[k][l]];
            _weights[k] *= w[n][indicesQ[k][l]];
        }
    }
}

std::vector<Vector> TensorizedQuadrature::GetNodes() {
    std::vector<Vector> nodes;
    for( unsigned i = 0; i < _nQTotal; ++i ) {
        Vector tmp( _numDimXi, 0.0 );
        for( unsigned j = 0; j < _numDimXi; ++j ) tmp[j] = _nodes[j][i];
        nodes.push_back( tmp );
    }
    return nodes;
}

Vector TensorizedQuadrature::GetWeights() {
    Vector weights( _nQTotal );
    for( unsigned i = 0; i < _nQTotal; ++i ) {
        weights[i] = _weights[i];
    }
    return weights;
}

unsigned TensorizedQuadrature::GetNodeCount() { return _nQTotal; }
