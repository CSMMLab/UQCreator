#include "sparsegrid.h"
//#include "normalsparsegrid.h"
#include "uniformsparsegrid.h"

SparseGrid::SparseGrid( unsigned dim, unsigned level ) : _dim( dim ), _level( level ), _nodeCount( 0 ), _unionCount( 0 ) {}

SparseGrid::~SparseGrid() {
    for( unsigned i = 0; i < _dim; i++ ) delete[] _nodes[i];
    delete[] _nodes;

    delete[] _weights;

    for( unsigned i = 0; i < _dim; i++ ) delete[] _unionIndex[i];
    delete[] _unionIndex;

    for( unsigned i = 0; i < _dim; i++ ) delete[] _nodeOrder[i];
    delete[] _nodeOrder;
}

///////////////////////////////////////////////////////////////////////////////////////

void SparseGrid::computeUnionIndex( unsigned dim, unsigned level, unsigned& nodeCount, unsigned& unionCount, unsigned**& unionIndex ) {
    std::vector<unsigned> level_mult( 1, 1 );
    for( unsigned i = 0; i < level; i++ ) {
        if( i != 1 ) {
            level_mult.push_back( level_mult.back() * 2 );
        }
        else {
            level_mult.push_back( level_mult.back() );
        }
    }

    std::vector<std::vector<unsigned>> uindex( dim );

    for( unsigned i = 0; i <= level; i++ ) {
        std::vector<unsigned> comp( dim, 0 );
        comp[0]      = i;
        unsigned tmp = comp[0];
        unsigned pos = 0;
        while( comp.back() != i ) {
            for( unsigned j = 0; j < dim; j++ ) {
                ( uindex[j] ).push_back( comp[j] );
            }
            if( tmp > 1 ) pos = 0;
            pos += 1;
            tmp           = comp[pos - 1];
            comp[pos - 1] = 0;
            comp[0]       = tmp - 1;
            comp[pos]     = comp[pos] + 1;

            std::vector<unsigned> prod;
            for( unsigned j = 0; j < dim; j++ ) {
                prod.push_back( level_mult[comp[j]] );
            }
            nodeCount += std::accumulate( prod.begin(), prod.end(), 1u, std::multiplies<unsigned>() );
        }
        for( unsigned j = 0; j < dim; j++ ) {
            ( uindex[j] ).push_back( comp[j] );
        }

        std::vector<unsigned> prod;
        for( unsigned j = 0; j < dim; j++ ) {
            prod.push_back( level_mult[comp[j]] );
        }
        nodeCount += std::accumulate( prod.begin(), prod.end(), 1u, std::multiplies<unsigned>() );
    }

    unionIndex = new unsigned*[dim];
    for( unsigned i = 0; i < dim; i++ ) {
        unionIndex[i] = new unsigned[( uindex[0] ).size()];
    }

    for( unsigned i = 0; i < dim; i++ )
        for( unsigned j = 0; j < uindex[0].size(); j++ ) unionIndex[i][j] = 0;

    for( unsigned i = 0; i < dim; i++ ) {
        for( unsigned j = 0; j < ( uindex[0] ).size(); j++ ) {
            unionIndex[i][j] = ( uindex[i] ).at( j );
        }
    }

    unionCount = static_cast<unsigned>( ( uindex[0] ).size() );
}

void SparseGrid::computeNodeOrder( unsigned dim,
                                   unsigned level,
                                   unsigned unionCount,
                                   unsigned** unionIndex,
                                   unsigned nodeCount,
                                   std::vector<std::vector<std::vector<unsigned>>>& gridIndices,
                                   unsigned**& nodeOrder ) {
    gridIndices =
        std::vector<std::vector<std::vector<unsigned>>>( unionCount, std::vector<std::vector<unsigned>>( dim, std::vector<unsigned>( 0 ) ) );
    std::vector<std::vector<unsigned>> grid_buffer( dim, std::vector<unsigned>( 0 ) );
    std::vector<unsigned> order_buffer( 0 );

    for( unsigned u = 0; u < unionCount; u++ ) {
        std::vector<unsigned> order1D( dim );

        for( unsigned i = 0; i < dim; i++ ) {
            if( unionIndex[i][u] == 0 )
                order1D[i] = 1;
            else
                order1D[i] = this->order( unionIndex[i][u] );
        }

        unsigned orderND = std::accumulate( order1D.begin(), order1D.end(), 1u, std::multiplies<unsigned>() );

        std::vector<std::vector<unsigned>> grid_index( dim, std::vector<unsigned>( orderND ) );
        std::vector<unsigned> a( dim, 0 );
        bool run_colex = true;
        unsigned p     = 0;
        while( true ) {
            colexOrdering( order1D, a, run_colex );
            if( !run_colex ) break;
            p++;
            for( unsigned i = 0; i < dim; i++ ) {
                grid_index[i][p] = a[i];
            }
        }

        for( unsigned i = 0; i < dim; i++ )
            for( unsigned j = 0; j < orderND; j++ ) gridIndices[u][i].push_back( grid_index[i][j] );

        unsigned order_max = 0;

        for( unsigned i = 0; i < dim; i++ ) {
            if( unionIndex[i][u] == 0 ) {
                if( level == 0 ) {
                    order_max = 1;
                }
                else {
                    order_max = this->order( level );
                }
                for( unsigned j = 0; j < orderND; j++ ) {
                    grid_index[i][j] = static_cast<unsigned>( std::floor( ( order_max - 1 ) / 2 ) );
                }
            }
            else {
                unsigned diff = this->order( level - unionIndex[i][u] ) - 1u;
                for( unsigned j = 0; j < orderND; j++ ) {
                    grid_index[i][j] *= diff;
                }
            }
        }

        for( unsigned i = 0; i < orderND; i++ ) {
            order_buffer.push_back( orderND );
            for( unsigned j = 0; j < dim; j++ ) grid_buffer[j].push_back( grid_index[j][i] );
        }
    }

    nodeOrder = new unsigned*[dim];
    for( unsigned i = 0; i < dim; i++ ) nodeOrder[i] = new unsigned[nodeCount];

    for( unsigned i = 0; i < dim; i++ )
        for( unsigned j = 0; j < nodeCount; j++ ) nodeOrder[i][j] = 0;

    std::vector<unsigned> maxNodeOrder( 0 );
    for( unsigned i = 0; i < nodeCount; i++ ) maxNodeOrder.push_back( std::numeric_limits<unsigned>::infinity() );

    unsigned ctr = 0;
    for( unsigned i = 0; i < grid_buffer[0].size(); i++ ) {
        bool found = false;
        if( ctr > nodeCount ) break;
        for( unsigned j = 0; j < ctr; j++ ) {
            bool match = true;
            for( unsigned k = 0; k < dim; k++ ) {
                if( nodeOrder[k][j] != grid_buffer[k][i] ) {
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
            for( unsigned k = 0; k < dim; k++ ) {
                nodeOrder[k][ctr] = grid_buffer[k][i];
            }
            ctr += 1;
        }
    }
}

template <class T> void SparseGrid::colexOrdering( std::vector<T> base, std::vector<T>& a, bool& run ) {
    for( unsigned i = 0; i < a.size(); ++i ) {
        a[i] += 1;
        if( a[i] < base[i] ) return;
        a[i] = 0;
    }
    run = false;
}

std::vector<double>
SparseGrid::prod( unsigned index, unsigned order, unsigned orderND, std::vector<double> weights1D, std::vector<double> weightsND ) {
    static unsigned contig;
    static double rep;
    static unsigned skip;

    if( index == 0 ) {
        contig = 1;
        skip   = 1;
        rep    = orderND;
        for( unsigned i = 0; i < weightsND.size(); i++ ) weightsND[i] = 1.0;
    }

    rep /= order;
    skip *= order;

    for( unsigned i = 0; i < order; i++ ) {
        unsigned start = i * contig;
        for( unsigned j = 0; j < rep; j++ ) {
            for( unsigned k = start; k < start + contig; k++ ) {
                weightsND[k] *= weights1D[i];
            }
            start = start + skip;
        }
    }

    contig *= order;

    return weightsND;
}

int SparseGrid::modu( unsigned i ) {
    if( i % 2 == 0 )
        return 1;
    else
        return -1;
}

int SparseGrid::choose( int n, int k ) {
    int min = std::min( k, n - k );
    if( min < 0 )
        return 0;
    else if( min == 0 )
        return 1;
    else {
        int max = std::max( k, n - k );
        int res = 1;
        for( int i = 1; i <= min; i++ ) res = ( res * ( max + i ) ) / i;
        return res;
    }
}

unsigned SparseGrid::GetNodeCount() { return _nodeCount; }

std::vector<Vector> SparseGrid::GetNodes() {
    std::cout << "Getting " << _nodeCount << " nodes" << std::endl;
    std::vector<Vector> nodes;
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        Vector tmp( _dim, 0.0 );
        for( unsigned j = 0; j < _dim; ++j ) tmp[j] = _nodes[j][i];
        nodes.push_back( tmp );
    }
    return nodes;
}

Vector SparseGrid::GetWeights() {
    Vector weights( _nodeCount );
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        weights[i] = _weights[i];
    }
    return weights;
}

std::vector<Vector> SparseGrid::GetNodeOrder() {
    std::vector<Vector> nodeOrder;
    for( unsigned i = 0; i < _nodeCount; ++i ) {
        Vector tmp( _dim, 0.0 );
        for( unsigned j = 0; j < _dim; ++j ) tmp[j] = _nodeOrder[j][i];
        nodeOrder.push_back( tmp );
    }
    return nodeOrder;
}

std::vector<Vector> SparseGrid::GetUnionIndex() {
    std::vector<Vector> unionIndex;
    for( unsigned i = 0; i < _unionCount; ++i ) {
        Vector tmp( _dim, 0.0 );
        for( unsigned j = 0; j < _dim; ++j ) tmp[j] = _unionIndex[j][i];
        unionIndex.push_back( tmp );
    }
    return unionIndex;
}

/*
SparseGrid* SparseGrid::Create( std::string type, unsigned dim, unsigned level ) {
    if( !( type.length() != 0 && dim > 0 ) ) {
        std::cerr << "\n[SparseGrid]: Invalid input parameter" << std::endl;
        exit( EXIT_FAILURE );
    }
    if( !std::strcmp( type.c_str(), "CC" ) || !std::strcmp( type.c_str(), "uniform" ) ) {
        return new UniformSparseGrid( dim, level );
    }
    else if( !std::strcmp( type.c_str(), "GH" ) || !std::strcmp( type.c_str(), "normal" ) ) {
        return new UniformSparseGrid( dim, level );
        // return new NormalSparseGrid( dim, level );
    }
    else {
        std::cerr << "\n[SparseGrid]: Invalid grid type" << std::endl;
        exit( EXIT_FAILURE );
    }
}*/
