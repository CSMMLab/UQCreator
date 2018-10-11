#include "mesh.h"

Mesh::Mesh( std::string inputFile ) : _status( MESH_STATUS_UNLOADED ) {
    auto file            = cpptoml::parse_file( inputFile );
    auto settings        = file->get_table( "mesh" );
    std::string meshFile = settings->get_as<std::string>( "fileName" ).value_or( "none" );
    if( meshFile.compare( "none" ) ) {
        this->Load( meshFile );
    }
    else {
        _dimension = settings->get_as<unsigned>( "Dimension" ).value_or( 0 );
        _numCells  = settings->get_as<unsigned>( "NumberOfCells" ).value_or( 0 );
        double a   = settings->get_as<double>( "a" ).value_or( 0.0 );
        double b   = settings->get_as<double>( "b" ).value_or( 0.0 );
        if( ( std::fabs( a ) <= std::numeric_limits<double>::epsilon() && std::fabs( b ) <= std::numeric_limits<double>::epsilon() ) ||
            _dimension == 0 || _numCells == 0 )
            std::cerr << "[ERROR]: Mesh::Mesh invalid mesh parameters" << std::endl;
        else {
            this->CreateGrid( a, b );
            _meshType = MESH_TYPE_1DPLAIN;
            _status   = MESH_STATUS_LOADED;
        }
        _neighbors.resize( _numCells + 4 );

        for( unsigned j = 2; j < _numCells + 1; ++j ) {

            _neighbors[j].resize( 2 );
            _neighbors[j][0] = j - 1;
            _neighbors[j][1] = j + 1;
        }
    }
}

Mesh::~Mesh() {}

void Mesh::Load( std::string filename ) { std::cerr << "[ERROR]: Mesh::Load not yet implemented" << std::endl; }

void Mesh::CreateGrid( double a, double b ) {
    assert( _dimension == 1 );
    _mesh.resize( _numCells + 4 );
    _spacing.resize( _numCells + 4 - 1 );

    for( unsigned j = 0; j < _numCells + 4; ++j ) {
        _mesh[j] = a + j * ( b - a ) / ( _numCells + 3 );
    }
    _n.resize( _numCells + 4 );
    _nUnit.resize( _numCells + 4 );
    Vector n1( 1 );
    Vector n2( 1 );
    for( unsigned j = 0; j <= _numCells + 2; ++j ) {
        _n[j].resize( 2 );
        _nUnit[j].resize( 2 );
        n1[0]        = -1.0;
        n2[0]        = 1.0;
        _nUnit[j][0] = n1;
        _nUnit[j][1] = n2;
        _spacing[j]  = _mesh[j + 1] - _mesh[j];
        _n[j][0]     = _nUnit[j][0] / _spacing[j];
        _n[j][1]     = _nUnit[j][1] / _spacing[j];
    }
}

unsigned Mesh::GetNumCells() const { return _numCells; }

unsigned Mesh::GetDimension() const { return _dimension; }

const Vector& Mesh::GetGrid() const { return _mesh; }

const Vector& Mesh::GetSpacing() const { return _spacing; }

double Mesh::GetArea( unsigned j ) const { return _spacing[j]; }

blaze::DynamicVector<unsigned> Mesh::GetNeighbors( unsigned j ) const { return _neighbors[j]; };

Vector Mesh::GetNormals( unsigned j, unsigned l ) const { return _n[j][l]; };

Vector Mesh::GetUnitNormals( unsigned j, unsigned l ) const { return _nUnit[j][l]; };
