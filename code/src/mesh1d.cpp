#include "mesh1d.h"

Mesh1D::Mesh1D( std::string inputFile ) : Mesh( 1 ) {
    auto file     = cpptoml::parse_file( inputFile );
    auto settings = file->get_table( "mesh" );
    _numCells     = settings->get_as<unsigned>( "NumberOfCells" ).value_or( 0 );
    double a      = settings->get_as<double>( "a" ).value_or( 0.0 );
    double b      = settings->get_as<double>( "b" ).value_or( 0.0 );
    if( ( std::fabs( a ) <= std::numeric_limits<double>::epsilon() && std::fabs( b ) <= std::numeric_limits<double>::epsilon() ) || _numCells == 0 )
        std::cerr << "[ERROR]: Mesh1D::Mesh invalid mesh parameters" << std::endl;
    else {
        this->CreateGrid( a, b );
    }
    _cells[0]->AddNeighbor( _cells[1] );
    for( unsigned i = 1; i < _numCells - 1; ++i ) {
        _cells[i]->AddNeighbor( _cells[i - 1] );
        _cells[i]->AddNeighbor( _cells[i + 1] );
    }
    _cells[_numCells - 1]->AddNeighbor( _cells[_numCells - 2] );

    _neighbors.resize( _numCells );
    _boundaryType.resize( _numCells );
    for( unsigned i = 1; i < _numCells - 1; ++i ) {
        _neighbors[i].resize( 2 );
        _neighbors[i][0] = i - 1;
        _neighbors[i][1] = i + 1;
        _boundaryType[i] = BoundaryType::NONE;
    }
    _neighbors[0].resize( 2 );
    _neighbors[0][0] = _numCells;
    _neighbors[0][1] = _numCells;
    _neighbors[_numCells - 1].resize( 2 );
    _neighbors[_numCells - 1][0] = _numCells;
    _neighbors[_numCells - 1][1] = _numCells;
    _boundaryType[0]             = BoundaryType::DIRICHLET;
    _boundaryType[_numCells - 1] = BoundaryType::DIRICHLET;
}

Mesh1D::~Mesh1D() {}

void Mesh1D::CreateGrid( double a, double b ) {
    _cells.resize( _numCells );
    _nodes.resize( _numCells + 1 );
    for( unsigned i = 0; i < _numCells + 1; ++i ) {
        if( i == 0 || i == _numCells ) {
            _nodes[i] = new Node{i, true, std::vector<double>( _dimension, a + i * ( b - a ) / ( _numCells ) )};
        }
        else {
            _nodes[i] = new Node{i, false, std::vector<double>( _dimension, a + i * ( b - a ) / ( _numCells ) )};
        }
    }
    for( unsigned i = 0; i < _numCells; ++i ) {
        std::vector<Node*> cellNodes{_nodes[i], _nodes[i + 1]};
        _cells[i] = new Line( i, cellNodes );
    }
}

void Mesh1D::Export() const {}

Vector Mesh1D::GetNodePositionsX() const {
    Vector x( _numCells, 0.0 );
    for( unsigned i = 0; i < _numCells; ++i ) {
        x[i] = _cells[i]->GetCenter()[0];
    }
    return x;
}
