#include "mesh2d.h"

Mesh2D::Mesh2D( Settings* settings ) : Mesh( settings, 2 ) {
    auto file         = cpptoml::parse_file( _settings->GetInputFile() );
    auto table        = file->get_table( "mesh" );
    auto formatString = table->get_as<std::string>( "format" );
    if( formatString ) {
        if( formatString->compare( "SU2" ) == 0 ) {
            _format = MeshFormat::SU2;
        }
        else {
            std::cerr << "[Inputfile][mesh][(string)format] Unsupported mesh format!\nPlease set one of the following formats: SU2" << std::endl;
            exit( EXIT_FAILURE );
        }
    }
    else {
        std::cerr << "[Inputfile][mesh][(string)format] Not set!\nPlease set one of the following formats: SU2" << std::endl;
        exit( EXIT_FAILURE );
    }

    if( _format == MeshFormat::SU2 ) {
        auto SU2MeshFile = table->get_as<std::string>( "SU2MeshFile" );
        if( SU2MeshFile ) {
            _SU2MeshFile = _settings->GetInputDir() + "/" + *SU2MeshFile;
        }
        else {
            std::cerr << "[Inputfile][mesh][(string)SU2MeshFile] Not set!" << std::endl;
            exit( EXIT_FAILURE );
        }
        auto BCStrings = table->get_array_of<cpptoml::array>( "SU2BC" );
        if( BCStrings ) {
            for( unsigned i = 0; i < BCStrings->size(); ++i ) {
                auto BCString = ( *BCStrings )[i]->get_array_of<std::string>();
                BoundaryType type;
                if( ( *BCString )[1].compare( "noslip" ) == 0 ) {
                    type = BoundaryType::NOSLIP;
                }
                else if( ( *BCString )[1].compare( "dirichlet" ) == 0 ) {
                    type = BoundaryType::DIRICHLET;
                }
                else if( ( *BCString )[1].compare( "neumann" ) == 0 ) {
                    type = BoundaryType::NEUMANN;
                }
                else if( ( *BCString )[1].compare( "periodic" ) == 0 ) {
                    type = BoundaryType::PERIODIC;
                }
                else {
                    std::cerr << "[Inputfile][mesh][SU2BC] Invalid boundary condition on boundary '" + ( *BCString )[0] + "'!" << std::endl;
                    exit( EXIT_FAILURE );
                }
                _BCs.push_back( std::make_pair( ( *BCString )[0], type ) );
            }
        }
        else {
            std::cerr << "[Inputfile][mesh][SU2BC] Not set!" << std::endl;
            exit( EXIT_FAILURE );
        }
        _outputFile = table->get_as<std::string>( "outputFile" ).value_or( "" );
        LoadSU2MeshFromFile( _SU2MeshFile );
    }
    _settings->SetNumCells( _numCells );
}

Mesh2D::~Mesh2D() {}

void Mesh2D::LoadSU2MeshFromFile( std::string meshfile ) {
    std::ifstream ifs( meshfile, std::ios::in );
    std::string line;
    if( ifs.is_open() ) {
        while( getline( ifs, line ) ) {
            if( line.find( "NDIME", 0 ) != std::string::npos ) {
                assert( _dimension == GetTrailingNumber( line ) );
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NPOIN", 0 ) != std::string::npos ) {
                unsigned numPoints = GetTrailingNumber( line );
                for( unsigned i = 0; i < numPoints; ++i ) {
                    getline( ifs, line );
                    std::stringstream ss;
                    ss << line;
                    unsigned id = 0;
                    double tmp;
                    std::vector<double> coords( 3, 0.0 );
                    while( !ss.eof() ) {
                        for( unsigned d = 0; d < _dimension; ++d ) {
                            ss >> tmp;
                            coords[d] = tmp;
                        }
                        ss >> id;
                    }
                    _nodes.push_back( new Node{id, false, coords} );
                }
                break;
            }
        }
        std::sort( _nodes.begin(), _nodes.end(), []( auto& left, auto& right ) { return left->id < right->id; } );

        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NMARK", 0 ) != std::string::npos ) {
                unsigned numBCs = GetTrailingNumber( line );
                for( unsigned i = 0; i < numBCs; ++i ) {
                    std::string markerTag;
                    std::vector<BoundaryElement> boundaryElements;
                    for( unsigned k = 0; k < 2; ++k ) {
                        getline( ifs, line );
                        if( line.find( "MARKER_TAG", 0 ) != std::string::npos ) {
                            markerTag    = line.substr( line.find( "=" ) + 1 );
                            auto end_pos = std::remove_if( markerTag.begin(), markerTag.end(), isspace );
                            markerTag.erase( end_pos, markerTag.end() );
                        }
                        else if( line.find( "MARKER_ELEMS", 0 ) != std::string::npos ) {
                            unsigned numMarkerElements = GetTrailingNumber( line );
                            for( unsigned j = 0; j < numMarkerElements; ++j ) {
                                getline( ifs, line );
                                std::stringstream ss;
                                ss << line;
                                unsigned type = 0, tmp = 0;
                                std::vector<unsigned> nodes;
                                while( !ss.eof() ) {
                                    ss >> type;
                                    for( unsigned d = 0; d < _dimension; ++d ) {
                                        ss >> tmp;
                                        nodes.push_back( tmp );
                                    }
                                }
                                for( auto& n : _nodes ) {
                                    for( auto& o : nodes ) {
                                        if( n->id == o ) {
                                            n->isBoundaryNode = true;
                                        }
                                    }
                                }
                                boundaryElements.push_back( BoundaryElement{type, nodes} );
                            }
                        }
                        else {
                            exit( EXIT_FAILURE );
                        }
                    }
                    BoundaryType type = BoundaryType::NONE;
                    for( const auto& b : _BCs ) {
                        if( b.first.compare( markerTag ) == 0 ) {
                            type = b.second;
                        }
                    }

                    assert( type != BoundaryType::NONE );
                    _boundaries.push_back( Boundary{markerTag, type, boundaryElements} );
                }
                break;
            }
        }
        ifs.clear();
        ifs.seekg( 0, std::ios::beg );
        while( getline( ifs, line ) ) {
            if( line.find( "NELEM", 0 ) != std::string::npos ) {
                unsigned numElements = GetTrailingNumber( line );
                _numCells            = numElements;
                for( unsigned i = 0; i < numElements; ++i ) {
                    getline( ifs, line );
                    std::stringstream ss;
                    ss << line;
                    unsigned type = 0, tmp = 0, id = 0;
                    std::vector<Node*> elementNodes;
                    while( !ss.eof() ) {
                        ss >> type;
                        unsigned nElementNodes;
                        if( type == 5 ) {
                            nElementNodes = 3;
                        }
                        else {
                            std::cerr << "Unsupported mesh type!" << std::endl;
                            exit( EXIT_FAILURE );
                        }
                        for( unsigned d = 0; d < nElementNodes; ++d ) {
                            ss >> tmp;
                            elementNodes.push_back( FindNodeByID( tmp ) );
                        }
                        ss >> id;
                    }
                    _cells.push_back( new Triangle( id, elementNodes ) );
                }
                break;
            }
        }
        _boundaryType.resize( _numCells );
        assert( _cells.size() == _numCells );
        for( unsigned i = 0; i < _numCells; ++i ) {
            if( _cells[i]->IsBoundaryCell() ) {
                for( unsigned k = 0; k < _cells[i]->GetNodeNum(); ++k ) {
                    unsigned nodeId = _cells[i]->GetNode( k )->id;
                    for( unsigned j = 0; j < _boundaries.size(); ++j ) {
                        for( unsigned l = 0; l < _boundaries[j].elements.size(); ++l ) {
                            if( _boundaries[j].elements[l].nodes[0] == nodeId || _boundaries[j].elements[l].nodes[1] == nodeId ) {
                                _boundaryType[i] = _boundaries[j].type;
                            }
                        }
                    }
                }
            }
        }
    }
    else {
        std::cerr << "File not found" << std::endl;
    }
    ifs.close();
    DetermineNeighbors();
    _neighborIDs.resize( _numCells );
    for( unsigned j = 0; j < _numCells; ++j ) {
        _neighborIDs[j] = _cells[j]->GetNeighborIDs();
        if( _cells[j]->IsBoundaryCell() ) {
            if( _neighborIDs[j][0] != _numCells && _neighborIDs[j][1] != _numCells && _neighborIDs[j][2] != _numCells ) {
                std::cerr << "Wrong boundary cell " << j << " detected" << std::endl;
                std::cerr << "Neighbors are " << _neighborIDs[j][0] << " " << _neighborIDs[j][1] << " " << _neighborIDs[j][2] << std::endl;
                std::cerr << "Points are " << _cells[j]->GetNodes()[0]->id << " " << _cells[j]->GetNodes()[1]->id << " "
                          << _cells[j]->GetNodes()[2]->id << std::endl;
                std::cerr << "Points A are " << _cells[_neighborIDs[j][0]]->GetNodes()[0]->id << " " << _cells[_neighborIDs[j][0]]->GetNodes()[1]->id
                          << " " << _cells[_neighborIDs[j][0]]->GetNodes()[2]->id << std::endl;
                std::cerr << "Points B are " << _cells[_neighborIDs[j][1]]->GetNodes()[0]->id << " " << _cells[_neighborIDs[j][1]]->GetNodes()[1]->id
                          << " " << _cells[_neighborIDs[j][1]]->GetNodes()[2]->id << std::endl;
                std::cerr << "Points C are " << _cells[_neighborIDs[j][2]]->GetNodes()[0]->id << " " << _cells[_neighborIDs[j][2]]->GetNodes()[1]->id
                          << " " << _cells[_neighborIDs[j][2]]->GetNodes()[2]->id << std::endl;
                exit( EXIT_FAILURE );
            }

            if( _cells[j]->GetNeighborIDs()[0] != _numCells && _cells[j]->GetNeighborIDs()[1] != _numCells &&
                _cells[j]->GetNeighborIDs()[2] != _numCells ) {
                std::cerr << "Wrong boundary cell " << j << " detected" << std::endl;
                exit( EXIT_FAILURE );
            }
        }
    }
}

unsigned Mesh2D::GetTrailingNumber( std::string s ) {
    return static_cast<unsigned>( std::stoi( s.substr( s.find_first_of( "0123456789" ), s.length() - 1 ) ) );
}

unsigned Mesh2D::BinarySearch( unsigned id, unsigned lBound, unsigned rBound ) {
    if( rBound >= lBound ) {
        unsigned mid = lBound + ( rBound - lBound ) / 2;
        if( _nodes[mid]->id == id ) return mid;
        if( _nodes[mid]->id > id ) return BinarySearch( id, lBound, mid - 1 );
        return BinarySearch( id, mid + 1, rBound );
    }
    std::cerr << "Binary search failed" << std::endl;
    exit( EXIT_FAILURE );
}

Node* Mesh2D::FindNodeByID( unsigned id ) {
    auto pos = BinarySearch( id, 0, static_cast<unsigned>( _nodes.size() - 1 ) );
    return _nodes[pos];
}

void Mesh2D::DetermineNeighbors() {
    unsigned index0 = 0, index1 = 0;
    for( auto& i : _cells ) {
        for( unsigned l = 0; l < i->GetNeighborIDs().size(); ++l ) {
            i->AddNeighborId( _numCells, 0 );
            i->AddNeighborId( _numCells, 1 );
            i->AddNeighborId( _numCells, 2 );
        }
        for( auto& j : _cells ) {
            // if( i->GetNeighbors().size() == i->GetNodeNum() - i->IsBoundaryCell() ) {
            //    goto cnt;
            //}
            if( i->GetID() != j->GetID() ) {
                unsigned matchCtr = 0;
                for( unsigned n = 0; n < i->GetNodeNum(); ++n ) {
                    Node* k = i->GetNodes()[n];
                    for( const auto& l : j->GetNodes() ) {
                        if( k->id == l->id ) {
                            matchCtr++;
                            if( matchCtr == 1 ) index0 = n;
                            if( matchCtr == 2 ) index1 = n;
                        }
                    }
                    if( matchCtr == 2 ) {
                        if( index0 == 0 && index1 == 2 )
                            i->AddNeighbor( j, 2 );
                        else
                            i->AddNeighbor( j, index0 );
                        break;
                    }
                }
            }
        }
        // cnt:;
    }
}

void Mesh2D::Export( Matrix results ) const {
    std::string vtkFile = _outputFile;
    auto writer         = vtkXMLUnstructuredGridWriterSP::New();
    if( vtkFile.substr( _outputFile.find_last_of( "." ) + 1 ) != "vtu" ) {
        vtkFile.append( "." );
        vtkFile.append( writer->GetDefaultFileExtension() );
    }
    writer->SetFileName( vtkFile.c_str() );
    auto grid = vtkUnstructuredGridSP::New();
    auto pts  = vtkPointsSP::New();
    pts->SetNumberOfPoints( static_cast<int>( _nodes.size() ) );
    for( const auto& node : _nodes ) {
        pts->SetPoint( node->id, node->coords[0], node->coords[1], node->coords[2] );
    }
    vtkCellArraySP cellArray = vtkCellArraySP::New();
    for( unsigned i = 0; i < _cells.size(); ++i ) {
        auto tri = vtkTriangleSP::New();
        for( unsigned j = 0; j < _cells[i]->GetNodeNum(); ++j ) {
            tri->GetPointIds()->SetId( j, _cells[i]->GetNode( j )->id );
        }
        cellArray->InsertNextCell( tri );
    }
    grid->SetCells( VTK_TRIANGLE, cellArray );

    auto cellData = vtkDoubleArraySP::New();
    cellData->SetName( "rho" );
    for( unsigned i = 0; i < _numCells; i++ ) {
        cellData->InsertNextValue( results( 0, i ) );
    }
    grid->GetCellData()->AddArray( cellData );

    cellData = vtkDoubleArraySP::New();
    cellData->SetName( "rhoU" );
    cellData->SetNumberOfComponents( 3 );
    cellData->SetComponentName( 0, "x" );
    cellData->SetComponentName( 1, "y" );
    cellData->SetComponentName( 2, "z" );
    cellData->SetNumberOfTuples( _numCells );
    for( unsigned i = 0; i < _numCells; i++ ) {
        cellData->SetTuple3( i, results( 1, i ), results( 2, i ), 0.0 );
    }
    grid->GetCellData()->AddArray( cellData );

    cellData = vtkDoubleArraySP::New();
    cellData->SetName( "rhoE" );
    for( unsigned i = 0; i < _numCells; i++ ) {
        cellData->InsertNextValue( results( 3, i ) );
    }
    grid->GetCellData()->AddArray( cellData );

    grid->SetPoints( pts );
    grid->Squeeze();

    auto converter = vtkCellDataToPointDataSP::New();
    converter->AddInputDataObject( grid );
    converter->PassCellDataOn();
    converter->Update();

    auto conv_grid = converter->GetOutput();

    writer->SetInputData( conv_grid );
    writer->SetDataModeToAscii();

    writer->Write();
}

Vector Mesh2D::GetNodePositionsX() const {
    Vector x( _numCells, 0.0 );
    for( unsigned i = 0; i < _numCells; ++i ) {
        x[i] = _cells[i]->GetCenter()[0];
    }
    return x;
}
