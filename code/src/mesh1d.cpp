#include "mesh1d.h"

#include <fstream>

Mesh1D::Mesh1D( Settings* settings ) : Mesh( settings, 1 ) {
    auto file     = cpptoml::parse_file( _settings->GetInputFile() );
    auto table    = file->get_table( "mesh" );
    auto numCells = table->get_as<unsigned>( "numberOfCells" );
    if( numCells ) {
        _numCells = *numCells;
        _settings->SetNumCells( _numCells );
    }
    else {
        _log->error( "[mesh1d] 'numberOfCells' not set!" );
    }
    double a = table->get_as<double>( "a" ).value_or( 0.0 );
    double b = table->get_as<double>( "b" ).value_or( 0.0 );
    if( ( std::fabs( a ) <= std::numeric_limits<double>::epsilon() && std::fabs( b ) <= std::numeric_limits<double>::epsilon() ) || _numCells == 0 )
        _log->error( "[mesh1d] Invalid mesh parameters!" );
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
    _cells[0]->SetBoundaryType( BoundaryType::DIRICHLET );
    _cells[_numCells - 1]->SetBoundaryType( BoundaryType::DIRICHLET );
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

std::vector<Vector> Mesh1D::Import() const {
    auto reader = vtkUnstructuredGridReaderSP::New();
    reader->SetFileName( _settings->GetICFile().c_str() );
    reader->Update();

    std::vector<Vector> data( _numCells, Vector( _settings->GetNStates() ) );
    auto grid     = reader->GetOutput();
    auto cellData = grid->GetCellData();
    for( unsigned i = 0; i < _numCells; ++i ) {
        if( _settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
            data[i][0] = cellData->GetArray( 0 )->GetTuple1( static_cast<int>( i ) );
        }
        else if( _settings->GetProblemType() == ProblemType::P_EULER_1D ) {
            data[i][0] = cellData->GetArray( 0 )->GetTuple1( static_cast<int>( i ) );
            data[i][1] = cellData->GetArray( 1 )->GetTuple3( static_cast<int>( i ) )[0];
            data[i][2] = cellData->GetArray( 2 )->GetTuple1( static_cast<int>( i ) );
        }
    }
    return data;
}

void Mesh1D::Export( const Matrix& results, std::string append ) const {
    assert( results.rows() >= _settings->GetNStates() * 2 );
    double height       = ( _nodes[_numCells]->coords[0] - _nodes[0]->coords[0] ) / 10;
    std::string vtkFile = _settings->GetOutputFile() + append;
    auto writer         = vtkUnstructuredGridWriterSP::New();
    if( vtkFile.substr( _outputFile.find_last_of( "." ) + 1 ) != "vtk" ) {
        vtkFile.append( ".vtk" );
    }
    writer->SetFileName( vtkFile.c_str() );
    auto grid = vtkUnstructuredGridSP::New();
    auto pts  = vtkPointsSP::New();
    pts->SetNumberOfPoints( static_cast<int>( _nodes.size() * 2 ) );
    for( const auto& node : _nodes ) {
        pts->SetPoint( node->id, node->coords[0], 0.0, 0.0 );
        pts->SetPoint( node->id + _numCells + 1, node->coords[0], height, 0.0 );
    }
    vtkCellArraySP cellArray = vtkCellArraySP::New();
    for( unsigned i = 0; i < _cells.size(); ++i ) {
        auto quad = vtkQuadSP::New();
        quad->GetPointIds()->SetId( 0, _cells[i]->GetNode( 0 )->id );
        quad->GetPointIds()->SetId( 1, _cells[i]->GetNode( 0 )->id + _numCells + 1 );
        quad->GetPointIds()->SetId( 2, _cells[i]->GetNode( 1 )->id + _numCells + 1 );
        quad->GetPointIds()->SetId( 3, _cells[i]->GetNode( 1 )->id );
        cellArray->InsertNextCell( quad );
    }
    grid->SetCells( VTK_QUAD, cellArray );

    if( _settings->GetProblemType() == ProblemType::P_BURGERS_1D ) {
        auto cellData = vtkDoubleArraySP::New();
        cellData->SetName( "E(u)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 0, i ) );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "Var(u)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 1, i ) );
        }
        grid->GetCellData()->AddArray( cellData );
    }
    else if( _settings->GetProblemType() == ProblemType::P_EULER_1D ) {
        auto cellData = vtkDoubleArraySP::New();
        cellData->SetName( "E(ρ)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 0, i ) );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "E(ρU)" );
        cellData->SetNumberOfComponents( 3 );
        cellData->SetComponentName( 0, "x" );
        cellData->SetComponentName( 1, "y" );
        cellData->SetComponentName( 2, "z" );
        cellData->SetNumberOfTuples( _numCells );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->SetTuple3( i, results( 1, i ), 0.0, 0.0 );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "E(ρE)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 2, i ) );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "Var(ρ)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 3, i ) );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "Var(ρU)" );
        cellData->SetNumberOfComponents( 3 );
        cellData->SetComponentName( 0, "x" );
        cellData->SetComponentName( 1, "y" );
        cellData->SetComponentName( 2, "z" );
        cellData->SetNumberOfTuples( _numCells );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->SetTuple3( i, results( 4, i ), 0.0, 0.0 );
        }
        grid->GetCellData()->AddArray( cellData );

        cellData = vtkDoubleArraySP::New();
        cellData->SetName( "Var(ρE)" );
        for( unsigned i = 0; i < _numCells; i++ ) {
            cellData->InsertNextValue( results( 5, i ) );
        }
        grid->GetCellData()->AddArray( cellData );
    }

    grid->SetPoints( pts );
    grid->Squeeze();

    auto converter = vtkCellDataToPointDataSP::New();
    converter->AddInputDataObject( grid );
    converter->PassCellDataOn();
    converter->Update();

    auto conv_grid = converter->GetOutput();

    writer->SetInputData( conv_grid );

    writer->Write();

    std::ofstream out( _settings->GetOutputFile() + "ExpectedValue" + append );
    for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
        out << GetCenterPos( j )[0];
        for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
            out << " " << results( s, j );
            if( _settings->HasExactSolution() ) out << " " << results( 2 * _settings->GetNStates() + s, j );
        }
        out << std::endl;
    }
    out.close();
    std::ofstream outV( _settings->GetOutputFile() + "Variance" + append );
    for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
        outV << GetCenterPos( j )[0];
        for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
            outV << " " << results( s + _settings->GetNStates(), j );
            if( results( s + _settings->GetNStates(), j ) < 0 ) {
                std::cerr << "Negative variance found in EXPORT" << std::endl;
                exit( EXIT_FAILURE );
            }
            if( _settings->HasExactSolution() ) outV << " " << results( 3 * _settings->GetNStates() + s, j );
        }
        outV << std::endl;
    }
    outV.close();
}

Vector Mesh1D::GetNodePositionsX() const {
    Vector x( _numCells, 0.0 );
    for( unsigned i = 0; i < _numCells; ++i ) {
        x[i] = _cells[i]->GetCenter()[0];
    }
    return x;
}
