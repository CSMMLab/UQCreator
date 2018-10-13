#include "mesh2d.h"

Mesh2D::Mesh2D( std::string inputFile ) : Mesh( 2 ) {
    auto file      = cpptoml::parse_file( inputFile );
    auto settings  = file->get_table( "mesh" );
    _SU2MeshFile   = settings->get_as<std::string>( "SU2File" ).value_or( "" );
    auto BCStrings = settings->get_array_of<cpptoml::array>( "SU2BC" );
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
            std::cerr << "Invalid boundary condition!" << std::endl;
        }
        _BCs.push_back( std::make_pair( ( *BCString )[0], type ) );
    }
    _outputFile = settings->get_as<std::string>( "outputFile" ).value_or( "" );
    LoadSU2MeshFromFile( _SU2MeshFile );
}

Mesh2D::~Mesh2D() {
    for( auto& i : _cells ) {
        delete i;
    }
    for( auto& i : _nodes ) {
        delete i;
    }
}

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
    }
    else {
        std::cerr << "File not found" << std::endl;
    }
    ifs.close();
    DetermineNeighbors();
}

void Mesh2D::ExportToVTK( std::string vtkfile ) const {
    auto writer = vtkXMLUnstructuredGridWriterSP::New();
    vtkfile.append( "." );
    vtkfile.append( writer->GetDefaultFileExtension() );
    writer->SetFileName( vtkfile.c_str() );
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
    grid->SetPoints( pts );
    grid->Squeeze();
    writer->SetInputData( grid );
    writer->SetDataModeToAscii();
    writer->Write();
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
    for( auto& i : _cells ) {
        for( auto& j : _cells ) {
            if( i->GetNeighbors().size() == i->GetNodeNum() - i->IsBoundaryCell() ) {
                goto cnt;
            }
            if( i->GetID() != j->GetID() ) {
                unsigned matchCtr = 0;
                for( const auto& k : i->GetNodes() ) {
                    for( const auto& l : j->GetNodes() ) {
                        if( k->id == l->id ) {
                            matchCtr++;
                        }
                    }
                    if( matchCtr == 2 ) {
                        i->AddNeighbor( j );
                        break;
                    }
                }
            }
        }
    cnt:;
    }
}

void Mesh2D::Export() const { ExportToVTK( _outputFile ); }
