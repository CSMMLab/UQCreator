#include "mesh1d.h"

Mesh1D::Mesh1D( Settings* settings ) : Mesh( settings, 1 ) {
    auto file     = cpptoml::parse_file( _settings->GetInputFile() );
    auto table    = file->get_table( "mesh" );
    auto numCells = table->get_as<unsigned>( "numberOfCells" );
    if( numCells ) {
        _numCells = *numCells;
        _settings->SetNumCells( _numCells );
    }
    else {
        std::cerr << "[Inputfile][mesh][(unsigned)numberOfCells] Not set!" << std::endl;
    }
    double a = table->get_as<double>( "a" ).value_or( 0.0 );
    double b = table->get_as<double>( "b" ).value_or( 0.0 );
    if( ( std::fabs( a ) <= std::numeric_limits<double>::epsilon() && std::fabs( b ) <= std::numeric_limits<double>::epsilon() ) || _numCells == 0 )
        std::cerr << "[Mesh1D]: Invalid mesh parameters!" << std::endl;
    else {
        this->CreateGrid( a, b );
    }
    _cells[0]->AddNeighbor( _cells[1], 1 );
    for( unsigned i = 1; i < _numCells - 1; ++i ) {
        _cells[i]->AddNeighbor( _cells[i - 1], 0 );
        _cells[i]->AddNeighbor( _cells[i + 1], 1 );
    }
    _cells[_numCells - 1]->AddNeighbor( _cells[_numCells - 2], 0 );

    _neighborIDs.resize( _numCells );
    _boundaryType.resize( _numCells );
    for( unsigned i = 1; i < _numCells - 1; ++i ) {
        _neighborIDs[i].resize( 2 );
        _neighborIDs[i][0] = i - 1;
        _neighborIDs[i][1] = i + 1;
        _boundaryType[i]   = BoundaryType::NONE;
    }
    _neighborIDs[0].resize( 2 );
    _neighborIDs[0][0] = _numCells;
    _neighborIDs[0][1] = _numCells;
    _neighborIDs[_numCells - 1].resize( 2 );
    _neighborIDs[_numCells - 1][0] = _numCells;
    _neighborIDs[_numCells - 1][1] = _numCells;
    _boundaryType[0]               = BoundaryType::DIRICHLET;
    _boundaryType[_numCells - 1]   = BoundaryType::DIRICHLET;
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

void Mesh1D::Export( Matrix results ) const {
    assert( results.rows() == _settings->GetNStates() * 2 );
    auto csvFile = _settings->GetOutputFile();
    if( csvFile.substr( _outputFile.find_last_of( "." ) + 1 ) != "csv" ) {
        csvFile.append( ".csv" );
    }
    std::ofstream writer( csvFile );
    for( unsigned i = 0; i < _settings->GetNStates() * 2; ++i ) {
        for( unsigned j = 0; j < _settings->GetNumCells() - 1; ++j ) {
            writer << results( i, j ) << ",";
        }
        writer << results( i, _settings->GetNumCells() - 1 ) << std::endl;
    }
    writer.close();
}

Vector Mesh1D::GetNodePositionsX() const {
    Vector x( _numCells, 0.0 );
    for( unsigned i = 0; i < _numCells; ++i ) {
        x[i] = _cells[i]->GetCenter()[0];
    }
    return x;
}
