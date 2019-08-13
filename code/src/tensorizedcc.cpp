#include "tensorizedcc.h"
#include "legendre.h"

TensorizedCC::TensorizedCC( Settings* settings ) : _settings( settings ), _nQuadPoints( settings->GetNQuadPoints() ) {
    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();
    if( _nQuadPoints == 0 )
        _nQTotal = 1;
    else
        _nQTotal = unsigned( std::pow( static_cast<unsigned>( std::pow( 2, _nQuadPoints ) + 1 ), _numDimXi ) );
    CreateGrid();
}

TensorizedCC::TensorizedCC( Settings* settings, unsigned nQuadPoints ) : _settings( settings ), _nQuadPoints( nQuadPoints ) {
    // compute total number of quad points
    _numDimXi = _settings->GetNDimXi();
    if( _nQuadPoints == 0 )
        _nQTotal = 1;
    else
        _nQTotal = unsigned( std::pow( static_cast<unsigned>( std::pow( 2, _nQuadPoints ) + 1 ), _numDimXi ) );
    CreateGrid();
}

TensorizedCC::~TensorizedCC() {
    for( unsigned i = 0; i < _numDimXi; i++ ) delete[] _nodes[i];
    delete[] _nodes;

    delete[] _weights;
}

void TensorizedCC::CreateGrid() {
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
    Vector xi = computeNodes1D( _nQuadPoints );
    Vector w  = computeWeights1D( _nQuadPoints );

    std::cout << "xi = " << xi << std::endl;
    std::cout << "w = " << w << std::endl;
    std::cout << std::endl;

    // setup weights and nodes
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            _nodes[l][k] = xi[indicesQ[k][l]];
            _weights[k] *= w[indicesQ[k][l]];
        }
    }
}

std::vector<Vector> TensorizedCC::GetNodes() {
    std::vector<Vector> nodes;
    for( unsigned i = 0; i < _nQTotal; ++i ) {
        Vector tmp( _numDimXi, 0.0 );
        for( unsigned j = 0; j < _numDimXi; ++j ) tmp[j] = _nodes[j][i];
        nodes.push_back( tmp );
    }
    return nodes;
}

Vector TensorizedCC::GetWeights() {
    Vector weights( _nQTotal );
    for( unsigned i = 0; i < _nQTotal; ++i ) {
        weights[i] = _weights[i];
    }
    return weights;
}

unsigned TensorizedCC::GetNodeCount() { return _nQTotal; }

/*void TensorizedCC::computeNodes( unsigned level, unsigned nodeCount, unsigned* nodeOrder, double*& nodes ) {
    unsigned n;
    unsigned Rl;
    if( level == 0 ) {
        Rl = 1;
    }
    else {
        Rl = static_cast<unsigned>( std::pow( 2, level ) + 1 );
    }
    for( unsigned j = 0; j < nodeCount; j++ ) {
        n = nodeOrder[j] + 1;
        if( n < 1 || Rl < n ) {
            nodes[j] = std::numeric_limits<unsigned>::min();
        }
        else if( Rl == 1 )
            nodes[j] = 0.0;
        else if( 2 * ( Rl - n ) == Rl - 1 )
            nodes[j] = 0.0;
        else
            nodes[j] = std::cos( ( Rl - n ) * M_PI / ( Rl - 1 ) );
    }
}*/

Vector TensorizedCC::computeNodes1D( unsigned level ) {
    double value;
    unsigned order;
    Vector nodes;
    if( level == 0 ) {
        order = 1;
    }
    else {
        order = static_cast<unsigned>( std::pow( 2, level ) + 1 );
    }
    nodes.resize( order );
    for( unsigned i = 1; i <= order; ++i ) {
        if( order < 1 )
            value = std::numeric_limits<unsigned>::min();
        else if( i < 1 || order < i )
            value = std::numeric_limits<unsigned>::min();
        else if( order == 1 )
            value = 0.0;
        else if( 2 * ( order - i ) == order - 1 )
            value = 0.0;
        else
            value = std::cos( ( order - i ) * M_PI / ( order - 1 ) );
        nodes[i - 1] = value;
    }

    return nodes;
}

Vector TensorizedCC::computeWeights1D( unsigned level ) {
    unsigned order;
    if( level == 0 ) {
        order = 1;
    }
    else {
        order = static_cast<unsigned>( std::pow( 2, level ) + 1 );
    }

    if( order == 1 ) {
        return Vector( 1, 2.0 );
    }

    Vector w( order, 0.0 );
    Vector theta( order, 0.0 );
    for( unsigned i = 0; i < order; i++ ) {
        theta[i] = i * M_PI / ( order - 1 );
    }

    for( unsigned i = 0; i < order; i++ ) {
        w[i] = 1.0;
        for( unsigned j = 1; j <= std::floor( ( order - 1 ) / 2 ); j++ ) {
            unsigned b;
            if( 2 * j == order - 1 )
                b = 1.0;
            else
                b = 2.0;
            w[i] -= b * std::cos( 2.0 * j * theta[i] ) / ( 4.0 * j * j - 1.0 );
        }
    }
    w[0] /= order - 1;
    for( unsigned i = 1; i < order - 1; i++ ) w[i] *= 2.0 / ( order - 1 );
    w[order - 1] /= order - 1;
    return w;
}
