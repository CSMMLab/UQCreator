#include "sparsegrid.h"

SparseGrid::SparseGrid( std::string type, unsigned pDim, unsigned pLevel ) : _dim( pDim ), _level( pLevel ), _nodeCount( 0 ), _unionCount( 0 ) {
    if( !( type.length() != 0 && pDim > 0 ) ) {
        std::cerr << "\n[sparseGrid]: Invalid input parameter" << std::endl;
        exit( EXIT_FAILURE );
    }
    if( _level == 0 ) {
        _Rl = 1;
    }
    else {
        _Rl = static_cast<unsigned>( std::pow( 2, _level ) + 1 );
    }

    if( !std::strcmp( type.c_str(), "CC" ) ) {
        createCCgrid();
    }
    else {
        std::cerr << "\n[sparseGrid]: Invalid node type" << std::endl;
        exit( EXIT_FAILURE );
    }
}

SparseGrid::~SparseGrid() {
    for( unsigned i = 0; i < _dim; i++ ) delete[] _nodes[i];
    delete[] _nodes;

    delete[] _weights;

    for( unsigned i = 0; i < _dim; i++ ) delete[] _unionIndex[i];
    delete[] _unionIndex;

    for( unsigned i = 0; i < _dim; i++ ) delete[] _nodeOrder[i];
    delete[] _nodeOrder;
}

void SparseGrid::createCCgrid() {
    computeUnionIndex();
    computeNodeOrder();

    _nodes = new double*[_dim];
    for( unsigned i = 0; i < _dim; i++ ) {
        _nodes[i] = new double[_nodeCount];
    }
    for( unsigned i = 0; i < _dim; i++ )
        for( unsigned j = 0; j < _nodeCount; j++ ) _nodes[i][j] = 0.0;
    computeCCNodes();

    _weights = new double[_nodeCount];
    for( unsigned i = 0; i < _nodeCount; i++ ) _weights[i] = 0.0;
    computeCCWeights();
}

void SparseGrid::computeCCNodes() {
    for( unsigned i = 0; i < _dim; i++ ) {
        for( unsigned j = 0; j < _nodeCount; j++ ) {
            unsigned n = _nodeOrder[i][j] + 1;
            if( n < 1 || _Rl < n ) {
                _nodes[i][j] = std::numeric_limits<int>::min();
            }
            else if( _Rl == 1 )
                _nodes[i][j] = 0.0;
            else if( 2 * ( _Rl - n ) == _Rl - 1 )
                _nodes[i][j] = 0.0;
            else
                _nodes[i][j] = std::cos( ( _Rl - n ) * M_PI / ( _Rl - 1 ) );
        }
    }
}

void SparseGrid::computeCCWeights() {
    if( _level == 0 ) {
        for( unsigned i = 0; i < _nodeCount; i++ ) _weights[i] = std::pow( 2, _dim );
        return;
    }

    unsigned min_level = std::max( 0u, _level + 1 - _dim );
    for( unsigned u = min_level; u < _unionCount; u++ ) {
        std::vector<unsigned> order1D( _dim );
        for( unsigned i = 0; i < _dim; i++ ) {
            if( _unionIndex[i][u] == 0 )
                order1D[i] = 1;
            else
                order1D[i] = static_cast<unsigned>( std::pow( 2, _unionIndex[i][u] ) + 1 );
        }

        unsigned orderND = std::accumulate( order1D.begin(), order1D.end(), 1u, std::multiplies<unsigned>() );

        std::vector<std::vector<unsigned>> grid_index = _gridIndices[u];

        std::vector<double> weightsND( orderND, 1.0 );
        for( unsigned i = 0; i < _dim; i++ ) {
            std::vector<double> weights1D = CCweights1D( order1D[i] );
            weightsND                     = CCweightsND( i, order1D[i], orderND, weights1D, weightsND );
        }

        unsigned l = 0;
        for( unsigned i = 0; i < _dim; i++ ) {
            if( _unionIndex[i][u] == 0 ) {
                unsigned order_max;
                if( _level == 0 )
                    order_max = 1;
                else
                    order_max = static_cast<unsigned>( std::pow( 2, _level ) + 1 );
                for( unsigned j = 0; j < orderND; j++ ) grid_index[i][j] = static_cast<unsigned>( std::floor( ( order_max - 1 ) / 2 ) );
            }
            else {
                unsigned factor = static_cast<unsigned>( std::pow( 2, _level - _unionIndex[i][u] ) );
                for( unsigned j = 0; j < orderND; j++ ) grid_index[i][j] *= factor;
            }
            l += _unionIndex[i][u];
        }
        if( l >= min_level ) {
            int coeff = static_cast<int>( isEven( _level - l ) * choose( _dim - 1u, _level - l ) );
#pragma omp parallel for schedule( dynamic, 2 )
            for( unsigned i = 0; i < orderND; i++ ) {
                for( unsigned j = 0; j < _nodeCount; j++ ) {
                    bool found = true;
                    for( unsigned k = 0; k < _dim; k++ ) {
                        if( grid_index[k][i] != _nodeOrder[k][j] ) {
                            found = false;
                            break;
                        }
                    }
                    if( found ) {
                        _weights[j] += coeff * weightsND[i];
                        break;
                    }
                }
            }
        }
    }
}

void SparseGrid::computeUnionIndex() {
    std::vector<unsigned> level_mult( 0 );
    level_mult.push_back( 1 );
    for( unsigned i = 0; i < _level; i++ ) {
        if( i != 1 ) {
            level_mult.push_back( level_mult.back() * 2 );
        }
        else {
            level_mult.push_back( level_mult.back() );
        }
    }

    std::vector<std::vector<unsigned>> uindex( _dim );

    for( unsigned i = 0; i <= _level; i++ ) {
        std::vector<unsigned> comp( _dim, 0 );
        comp[0]      = i;
        unsigned tmp = comp[0];
        unsigned pos = 0;
        while( comp.back() != i ) {
            for( unsigned j = 0; j < _dim; j++ ) {
                ( uindex[j] ).push_back( comp[j] );
            }
            if( tmp > 1 ) pos = 0;
            pos += 1;
            tmp           = comp[pos - 1];
            comp[pos - 1] = 0;
            comp[0]       = tmp - 1;
            comp[pos]     = comp[pos] + 1;

            std::vector<unsigned> prod;
            for( unsigned j = 0; j < _dim; j++ ) {
                prod.push_back( level_mult[comp[j]] );
            }
            _nodeCount += std::accumulate( prod.begin(), prod.end(), 1u, std::multiplies<unsigned>() );
        }
        for( unsigned j = 0; j < _dim; j++ ) {
            ( uindex[j] ).push_back( comp[j] );
        }

        std::vector<unsigned> prod;
        for( unsigned j = 0; j < _dim; j++ ) {
            prod.push_back( level_mult[comp[j]] );
        }
        _nodeCount += std::accumulate( prod.begin(), prod.end(), 1u, std::multiplies<unsigned>() );
    }

    _unionIndex = new unsigned*[_dim];
    for( unsigned i = 0; i < _dim; i++ ) {
        _unionIndex[i] = new unsigned[( uindex[0] ).size()];
    }

    for( unsigned i = 0; i < _dim; i++ )
        for( unsigned j = 0; j < uindex[0].size(); j++ ) _unionIndex[i][j] = 0;

    for( unsigned i = 0; i < _dim; i++ ) {
        for( unsigned int j = 0; j < ( uindex[0] ).size(); j++ ) {
            _unionIndex[i][j] = ( uindex[i] ).at( j );
        }
    }

    _unionCount = ( uindex[0] ).size();
}

void SparseGrid::computeNodeOrder() {
    _gridIndices =
        std::vector<std::vector<std::vector<unsigned>>>( _unionCount, std::vector<std::vector<unsigned>>( _dim, std::vector<unsigned>( 0 ) ) );
    std::vector<std::vector<unsigned>> grid_buffer( _dim, std::vector<unsigned>( 0 ) );
    std::vector<unsigned> order_buffer( 0 );

    for( unsigned u = 0; u < _unionCount; u++ ) {
        std::vector<unsigned> order1D( _dim );

        for( unsigned i = 0; i < _dim; i++ ) {
            if( _unionIndex[i][u] == 0 )
                order1D[i] = 1;
            else
                order1D[i] = static_cast<unsigned>( std::pow( 2, _unionIndex[i][u] ) + 1 );
        }

        unsigned orderND = std::accumulate( order1D.begin(), order1D.end(), 1u, std::multiplies<unsigned>() );

        std::vector<std::vector<unsigned>> grid_index( _dim, std::vector<unsigned>( orderND ) );
        std::vector<unsigned> a( _dim, 0 );
        bool run_colex = true;
        unsigned p     = 0;
        while( true ) {
            colexOrdering( order1D, a, run_colex );
            if( !run_colex ) break;
            p++;
            for( unsigned i = 0; i < _dim; i++ ) {
                grid_index[i][p] = a[i];
            }
        }

        for( unsigned i = 0; i < _dim; i++ )
            for( unsigned j = 0; j < orderND; j++ ) _gridIndices[u][i].push_back( grid_index[i][j] );

        unsigned order_max = 0;

        for( unsigned i = 0; i < _dim; i++ ) {
            if( _unionIndex[i][u] == 0 ) {
                if( _level == 0 ) {
                    order_max = 1;
                }
                else {
                    order_max = static_cast<unsigned>( std::pow( 2, _level ) + 1 );
                }
                for( unsigned j = 0; j < orderND; j++ ) {
                    grid_index[i][j] = static_cast<unsigned>( std::floor( ( order_max - 1 ) / 2 ) );
                }
            }
            else {
                unsigned diff = static_cast<unsigned>( std::pow( 2, _level - _unionIndex[i][u] ) );
                for( unsigned j = 0; j < orderND; j++ ) {
                    grid_index[i][j] *= diff;
                }
            }
        }

        for( unsigned i = 0; i < orderND; i++ ) {
            order_buffer.push_back( orderND );
            for( unsigned j = 0; j < _dim; j++ ) grid_buffer[j].push_back( grid_index[j][i] );
        }
    }

    _nodeOrder = new unsigned*[_dim];
    for( unsigned i = 0; i < _dim; i++ ) _nodeOrder[i] = new unsigned[_nodeCount];

    for( unsigned i = 0; i < _dim; i++ )
        for( unsigned j = 0; j < _nodeCount; j++ ) _nodeOrder[i][j] = 0;

    std::vector<unsigned> maxNodeOrder( 0 );
    for( unsigned i = 0; i < _nodeCount; i++ ) maxNodeOrder.push_back( std::numeric_limits<int>::infinity() );

    unsigned ctr = 0;
    for( unsigned i = 0; i < grid_buffer[0].size(); i++ ) {
        bool found = false;
        if( ctr > _nodeCount ) break;
        for( unsigned j = 0; j < ctr; j++ ) {
            bool match = true;
            for( unsigned k = 0; k < _dim; k++ ) {
                if( _nodeOrder[k][j] != grid_buffer[k][i] ) {
                    match = false;
                    break;
                }
            }
            if( match ) {
                found = true;
                if( maxNodeOrder[j] > order_buffer[i] ) {
                    maxNodeOrder[j] = order_buffer[i];
                }
                break;
            }
        }
        if( !found ) {
            for( unsigned k = 0; k < _dim; k++ ) {
                _nodeOrder[k][ctr] = grid_buffer[k][i];
            }
            ctr += 1;
        }
    }
}

std::vector<std::pair<std::vector<double>, double>> SparseGrid::GetGrid() const {
    std::vector<std::pair<std::vector<double>, double>> grid;
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        std::vector<double> tmp( _dim, 0.0 );
        for( unsigned j = 0; j < _dim; ++j ) tmp[j] = _nodes[j][i];
        grid.push_back( std::make_pair( tmp, _weights[i] ) );
    }
    return grid;
}

std::vector<std::vector<double>> SparseGrid::GetNodes() const {
    std::vector<std::vector<double>> nodes;
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        std::vector<double> tmp( _dim, 0.0 );
        for( unsigned j = 0; j < _dim; ++j ) tmp[j] = _nodes[j][i];
        nodes.push_back( tmp );
    }
    return nodes;
}
std::vector<double> SparseGrid::GetWeights() const {
    std::vector<double> weights;
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        weights.push_back( _weights[i] );
    }
    return weights;
}

unsigned long SparseGrid::GetNodeCount() const { return _nodeCount; }

void SparseGrid::colexOrdering( std::vector<unsigned> base, std::vector<unsigned>& a, bool& run ) {
    for( unsigned i = 0; i < _dim; i++ ) {
        a[i] += 1;
        if( a[i] < base[i] ) return;
        a[i] = 0;
    }
    run = false;
}

std::vector<double> SparseGrid::CCweights1D( unsigned order ) {
    if( order == 1 ) {
        return std::vector<double>( 1, 2.0 );
    }

    std::vector<double> w( order, 0.0 );
    std::vector<double> theta( order, 0.0 );
    for( unsigned i = 0; i < order; i++ ) {
        theta[i] = i * M_PI / ( order - 1 );
    }

    for( unsigned i = 0; i < order; i++ ) {
        w[i] = 1.0;
        for( unsigned j = 1; j <= std::floor( ( order - 1 ) / 2 ); j++ ) {
            int b;
            if( 2 * j == ( order - 1 ) )
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

std::vector<double>
SparseGrid::CCweightsND( unsigned dimIndex, unsigned order, unsigned orderND, std::vector<double> weights1D, std::vector<double> weightsND ) {
    static unsigned contig;
    static double rep;
    static unsigned skip;

    if( dimIndex == 0 ) {
        contig = 1;
        skip   = 1;
        rep    = orderND;
        for( unsigned i = 0; i < weightsND.size(); i++ ) weightsND[i] = 1.0;
    }

    rep /= order;
    skip *= order;

    for( unsigned i = 0; i < order; i++ ) {
        unsigned start = i * contig;
        for( int j = 0; j < rep; j++ ) {
            for( unsigned k = start; k < start + contig; k++ ) {
                weightsND[k] *= weights1D[i];
            }
            start = start + skip;
        }
    }
    contig *= order;
    return weightsND;
}

template <class T> int SparseGrid::isEven( T i ) {
    if( i % 2 == 0 )
        return 1;
    else
        return -1;
}

template <class T> double SparseGrid::choose( T n, T k ) {
    assert( n >= k );
    T min = std::min( k, n - k );
    if( min < 0 )
        return 0;
    else if( min == 0 )
        return 1;
    else {
        T max      = std::max( k, n - k );
        double res = 1;
        for( T i = 1; i <= min; i++ ) res = ( res * ( max + i ) ) / i;
        return res;
    }
}
