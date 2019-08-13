#include "tensorizedcc.h"
#include "legendre.h"
#include "uniformsparsegrid.h"

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

    unsigned nQuadLevel = static_cast<unsigned>( std::pow( 2, _nQuadPoints ) + 1 );
    unsigned currentIndex;
    std::vector<unsigned> tmp;
    tmp.resize( _numDimXi );
    bool push;

    // setup map from k (0,...,_nQTotal-1) to individual indices
    std::vector<std::vector<unsigned>> indicesQ;
    for( unsigned lev = 1; lev <= _nQuadPoints; ++lev ) {
        for( unsigned k = 0; k < _nQTotal; ++k ) {
            push = true;
            for( unsigned l = 0; l < _numDimXi; ++l ) {
                tmp[l] = unsigned( ( k - k % unsigned( std::pow( nQuadLevel, l ) ) ) / unsigned( std::pow( nQuadLevel, l ) ) ) % nQuadLevel;
                if( tmp[l] > pow( 2, lev ) ) {
                    push = false;
                    break;
                }
            }
            if( lev > 1 ) {
                unsigned max = 0;
                for( unsigned l = 0; l < _numDimXi; ++l ) {
                    if( max < tmp[l] ) max = tmp[l];
                }
                if( max < pow( 2, lev - 1 ) + 1 ) {
                    push = false;
                }
            }

            if( push ) indicesQ.push_back( tmp );
        }
    }

    auto quadGrid = new UniformSparseGrid( _nQuadPoints, 1 );

    // store quad points and weights for individual distributions
    std::vector<Vector> xi = quadGrid->GetNodes();
    std::cout << "size = " << xi.size() << std::endl;
    Vector w = quadGrid->GetWeights();

    // std::cout << "xi = " << xi << std::endl;
    // std::cout << "w = " << w << std::endl;
    // std::cout << std::endl;

    // setup weights and nodes
    for( unsigned k = 0; k < _nQTotal; ++k ) {
        for( unsigned l = 0; l < _numDimXi; ++l ) {
            _nodes[l][k] = xi[indicesQ[k][l]][0];
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
    unsigned counter = 3;
    _counter.resize( level );
    Vector nodes;
    if( level == 0 ) {
        order = 1;
        return Vector( 1, 0.0 );
    }
    else {
        order = static_cast<unsigned>( std::pow( 2, level ) + 1 );
    }
    nodes.resize( order );
    for( unsigned l = 1; l <= level; ++l ) {
        // std::cout << "Level " << l << std::endl;
        // determine order, i.e. number of points
        order = static_cast<unsigned>( std::pow( 2, l ) + 1 );

        if( l == 1 ) {
            nodes[0] = -1.0;
            _counter[l - 1].push_back( 0 );
            nodes[1] = 0.0;
            _counter[l - 1].push_back( 1 );
            nodes[2] = 1.0;
            _counter[l - 1].push_back( 2 );
            continue;
        }

        for( unsigned i = 2; i <= order; i += 2 ) {
            // std::cout << "Order " << i << std::endl;
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
            nodes[counter] = value;
            //_counter[l - 1][i - 1] = counter;
            ++counter;
        }
    }
    /*
        for( unsigned i = 0; i < _counter.size(); ++i ) {
            for( unsigned k = 0; k < _counter[i].size(); ++k ) {
                std::cout << _counter[i][k] << " ";
            }
            std::cout << std::endl;
        }*/

    return nodes;
}
/*
void TensorizedCC::DetermineCounter( unsigned level ) {
    unsigned order;
    unsigned counter = 3;
    _counter.resize( level );
    if( level == 0 ) {
        _counter[0].push_back( 0 );
        return;
    }

    for( unsigned l = 1; l <= level; ++l ) {
        order = static_cast<unsigned>( std::pow( 2, l ) + 1 );
        _counter[l - 1].resize( order );

        if( l == 1 ) {
            nodes[0] = -1.0;
            nodes[1] = 0.0;
            nodes[2] = 1.0;
            _counter[l - 1][];
            continue;
        }

        for( unsigned i = 2; i <= order; i += 2 ) {
            _counter[l - 1][i - 1] = counter;
            ++counter;
        }
    }

    return nodes;
}*/

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

    Vector wNew( order, 0.0 );

    for( unsigned l = 2; l <= level; ++l ) {
        for( unsigned j = pow( 2, level - l ) - 1; j <= order; j += std::floor( ( order - 1 ) / 2 * ( l - 1 ) ) ) {
            std::cout << j << " ";
        }
        std::cout << std::endl;
    }

    return w;
}
