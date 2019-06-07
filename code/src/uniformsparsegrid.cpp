#include "uniformsparsegrid.h"

UniformSparseGrid::UniformSparseGrid( Settings* settings ) : SparseGrid( settings->GetNDimXi(), settings->GetNQuadPoints() ) { createGrid(); }

UniformSparseGrid::~UniformSparseGrid() {}

///////////////////////////////////////////////////////////////////////////////////////

unsigned UniformSparseGrid::order( unsigned level ) { return static_cast<unsigned>( std::pow( 2, level ) + 1u ); }

void UniformSparseGrid::createGrid() {
    computeUnionIndex( _dim, _level, _nodeCount, _unionCount, _unionIndex );
    computeNodeOrder( _dim, _level, _unionCount, _unionIndex, _nodeCount, _gridIndices, _nodeOrder );

    _nodes = new double*[_dim];
    for( unsigned i = 0; i < _dim; i++ ) {
        _nodes[i] = new double[_nodeCount];
    }
    for( unsigned i = 0; i < _dim; i++ )
        for( unsigned j = 0; j < _nodeCount; j++ ) _nodes[i][j] = 0.0;
    computeNodes( _dim, _level, _nodeCount, _nodeOrder, _nodes );

    _weights = new double[_nodeCount];
    for( unsigned i = 0; i < _nodeCount; i++ ) _weights[i] = 0.0;
    computeWeights( _dim, _level, _nodeCount, _unionCount, _nodeOrder, _unionIndex, _weights );
}

void UniformSparseGrid::computeNodes( unsigned dim, unsigned level, unsigned nodeCount, unsigned** nodeOrder, double**& nodes ) {
    unsigned n;
    unsigned Rl;
    if( level == 0 ) {
        Rl = 1;
    }
    else {
        Rl = static_cast<unsigned>( std::pow( 2, level ) + 1 );
    }
    for( unsigned i = 0; i < dim; i++ ) {
        for( unsigned j = 0; j < nodeCount; j++ ) {
            n = nodeOrder[i][j] + 1;
            if( n < 1 || Rl < n ) {
                nodes[i][j] = std::numeric_limits<unsigned>::min();
            }
            else if( Rl == 1 )
                nodes[i][j] = 0.0;
            else if( 2 * ( Rl - n ) == Rl - 1 )
                nodes[i][j] = 0.0;
            else
                nodes[i][j] = std::cos( ( Rl - n ) * M_PI / ( Rl - 1 ) );
        }
    }
}

void UniformSparseGrid::computeWeights(
    unsigned dim, unsigned level, unsigned nodeCount, unsigned unionCount, unsigned** nodeOrder, unsigned** unionIndex, double*& weights ) {
    if( level == 0 ) {
        for( unsigned i = 0; i < nodeCount; i++ ) weights[i] = std::pow( 2, dim );
        return;
    }

    unsigned min_level = static_cast<unsigned>( std::max( 0, static_cast<int>( level ) + 1 - static_cast<int>( dim ) ) );
    for( unsigned u = min_level; u < unionCount; u++ ) {
        std::vector<unsigned> order1D( dim );
        for( unsigned i = 0; i < dim; i++ ) {
            if( unionIndex[i][u] == 0 )
                order1D[i] = 1;
            else
                order1D[i] = static_cast<unsigned>( std::pow( 2, unionIndex[i][u] ) + 1 );
        }

        unsigned orderND = std::accumulate( order1D.begin(), order1D.end(), 1u, std::multiplies<unsigned>() );

        std::vector<std::vector<unsigned>> grid_index = _gridIndices[u];

        std::vector<double> weightsND( orderND, 1.0 );
        for( unsigned i = 0; i < dim; i++ ) {
            std::vector<double> weights1D = computeWeights1D( order1D[i] );
            weightsND                     = prod( i, order1D[i], orderND, weights1D, weightsND );
        }

        unsigned l = 0;
        for( unsigned i = 0; i < dim; i++ ) {
            if( unionIndex[i][u] == 0 ) {
                unsigned order_max;
                if( level == 0 )
                    order_max = 1;
                else
                    order_max = static_cast<unsigned>( std::pow( 2, level ) + 1 );
                for( unsigned j = 0; j < orderND; j++ ) grid_index[i][j] = static_cast<unsigned>( std::floor( ( order_max - 1 ) / 2 ) );
            }
            else {
                unsigned factor = static_cast<unsigned>( std::pow( 2, level - unionIndex[i][u] ) );
                for( unsigned j = 0; j < orderND; j++ ) grid_index[i][j] *= factor;
            }
            l += unionIndex[i][u];
        }
        if( l >= min_level ) {
            int coeff = modu( level - l ) * choose( static_cast<int>( dim ) - 1, static_cast<int>( level ) - static_cast<int>( l ) );
            std::cout << modu( level - l ) << std::endl;
            for( unsigned i = 0; i < orderND; i++ ) {
                for( unsigned j = 0; j < nodeCount; j++ ) {
                    bool found = true;
                    for( unsigned k = 0; k < dim; k++ ) {
                        if( grid_index[k][i] != nodeOrder[k][j] ) {
                            found = false;
                            break;
                        }
                    }
                    if( found ) {
                        weights[j] += coeff * weightsND[i];
                        break;
                    }
                }
            }
        }
    }
    // for( unsigned i = 0; i < nodeCount; ++i ) weights[i] /= std::pow( 2, dim );
}

std::vector<double> UniformSparseGrid::computeWeights1D( unsigned order ) {
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
