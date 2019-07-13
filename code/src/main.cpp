#include <chrono>
#include <iostream>
#include <mpi.h>
#include <omp.h>

#include "mesh.h"
#include "momentsolver.h"
#include "problem.h"
#include "quadraturegrid.h"
#include "settings.h"
#include "tensorizedquadrature.h"
#include "uniformsparsegrid.h"

#include "spdlog/spdlog.h"

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"

Vector CalculateError( Settings* settings,
                       Mesh* mesh,
                       std::vector<std::vector<double>> referenceSolution,
                       const Matrix& solution,
                       unsigned LNorm,
                       const Vector& a,
                       const Vector& b ) {
    unsigned nStates = settings->GetNStates();
    unsigned nCells  = settings->GetNumCells();
    Vector error( 2 * nStates, 0.0 );
    Vector refNorm( 2 * nStates, 0.0 );
    if( settings->GetMeshDimension() == 2 ) {
        for( unsigned j = 0; j < nCells; ++j ) {
            if( mesh->GetGrid()[j]->GetCenter()[0] > a[0] && mesh->GetGrid()[j]->GetCenter()[0] < b[0] && mesh->GetGrid()[j]->GetCenter()[1] > a[1] &&
                mesh->GetGrid()[j]->GetCenter()[1] < b[1] ) {

                switch( LNorm ) {
                    case 1:
                        for( unsigned s = 0; s < 2 * nStates; ++s ) {
                            error[s] += std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j );
                            refNorm[s] += std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j );
                        }
                        break;
                    case 2:
                        for( unsigned s = 0; s < 2 * nStates; ++s ) {
                            error[s] += std::pow( ( solution( s, j ) - referenceSolution[j][s] ), 2 ) * mesh->GetArea( j );
                            refNorm[s] += std::pow( referenceSolution[j][s], 2 ) * mesh->GetArea( j );
                        }

                        break;
                    default: exit( EXIT_FAILURE );
                }
            }
        }
    }
    else {
        for( unsigned j = 0; j < nCells; ++j ) {
            switch( LNorm ) {
                case 1:
                    for( unsigned s = 0; s < 2 * nStates; ++s ) {
                        error[s] += std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j );
                        refNorm[s] += std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j );
                    }
                    break;
                case 2:
                    for( unsigned s = 0; s < 2 * nStates; ++s ) {
                        error[s] += std::pow( ( solution( s, j ) - referenceSolution[j][s] ), 2 ) * mesh->GetArea( j );
                        refNorm[s] += std::pow( referenceSolution[j][s], 2 ) * mesh->GetArea( j );
                    }

                    break;
                default: exit( EXIT_FAILURE );
            }
        }
    }
    for( unsigned s = 0; s < 2 * nStates; ++s ) {
        error[s] = std::pow( error[s] / refNorm[s], 1.0 / double( LNorm ) );
    }
    return error;
}

bool CheckInput( std::string& configFile, int argc, char* argv[] ) {
    std::string usage_help = "\n"
                             "Usage: " +
                             std::string( argv[0] ) +
                             " -c inputfile\n\n"
                             "Options:\n"
                             "  -t N             number of threads to use\n"
                             "  -h               displays this message\n";

    if( argc < 3 ) {
        std::cout << usage_help;
        return false;
    }
    for( int i = 1; i < argc; i++ ) {
        std::string arg = argv[i];
        if( arg == "-h" ) {
            std::cout << usage_help;
            return false;
        }
        else if( arg == "-c" ) {
            configFile = std::string( argv[++i] );
            std::ifstream f( configFile );
            if( !f.is_open() ) {
                std::cerr << "[ERROR] Unable to open specified inputfile!" << std::endl;
                return false;
            }
        }
        else if( arg == "-t" ) {
            omp_set_num_threads( std::stoi( argv[++i] ) );
        }
        else {
            std::cout << usage_help;
            return false;
        }
    }
    return true;
}

const std::string currentDateTime() {
    time_t now = time( nullptr );
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime( &now );
    strftime( buf, sizeof( buf ), "%Y-%m-%d_%X", &tstruct );
    return buf;
}

void initLogger( spdlog::level::level_enum terminalLogLvl, spdlog::level::level_enum fileLogLvl, std::string configFile ) {
    // event logger
    auto terminalSink = std::make_shared<spdlog::sinks::stdout_sink_mt>();
    terminalSink->set_level( terminalLogLvl );
    terminalSink->set_pattern( "%v" );

    auto file      = cpptoml::parse_file( configFile );
    auto general   = file->get_table( "general" );
    auto outputDir = general->get_as<std::string>( "outputDir" ).value_or( "." );
    if( !std::filesystem::exists( outputDir ) ) {
        std::filesystem::create_directory( outputDir );
    }
    if( !std::filesystem::exists( outputDir + "/logs" ) ) {
        std::filesystem::create_directory( outputDir + "/logs" );
    }

    int mype;
    MPI_Comm_rank( MPI_COMM_WORLD, &mype );
    std::ostringstream strs;
    strs << mype;
    std::string strPE = strs.str();

    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + currentDateTime() + "_PE" + strPE );
    fileSink->set_level( fileLogLvl );
    fileSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );

    std::vector<spdlog::sink_ptr> sinks;
    if( mype == 0 ) sinks.push_back( terminalSink );
    if( fileLogLvl != spdlog::level::off ) sinks.push_back( fileSink );

    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );

    auto momentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + currentDateTime() + "_moments" );
    momentFileSink->set_level( spdlog::level::info );
    momentFileSink->set_pattern( "%v" );
    auto moment_logger = std::make_shared<spdlog::logger>( "moments", momentFileSink );
    spdlog::register_logger( moment_logger );
}

void PrintInit( std::string configFile ) {
    auto log = spdlog::get( "event" );
    log->info( "UQCreator" );
    log->info( "==================================" );
    log->info( "" );
    log->info( "Config file:\t{0}", configFile );
    log->info( "==================================" );
    std::ifstream ifs( configFile );
    if( ifs.is_open() ) {
        std::string line;
        while( !ifs.eof() ) {
            std::getline( ifs, line );
            log->info( " {0}", line );
        }
    }
    log->info( "==================================" );
    log->info( "" );
}

int main( int argc, char* argv[] ) {
    MPI_Init( &argc, &argv );

    std::string configFile = "";

    if( !CheckInput( configFile, argc, argv ) ) {
        return EXIT_FAILURE;
    }

    initLogger( spdlog::level::info, spdlog::level::info, configFile );
    auto log = spdlog::get( "event" );

    // PrintInit( configFile );

    Settings* settings   = new Settings( configFile );
    QuadratureGrid* quad = QuadratureGrid::Create( settings );
    settings->SetNQTotal( quad->GetNodeCount() );
    if( settings->GetMyPE() == 0 ) PrintInit( configFile );
    Mesh* mesh           = Mesh::Create( settings );
    Problem* problem     = Problem::Create( settings );
    MomentSolver* solver = new MomentSolver( settings, mesh, problem );
    Closure* closure     = Closure::Create( settings );

    std::vector<MatVec> uQ;
    uQ.resize( settings->GetKEnd() - settings->GetKStart() + 1 );
    std::vector<Vector> xi = quad->GetNodes();
    unsigned nCells        = settings->GetNumCells();
    unsigned nStates       = settings->GetNStates();
    unsigned nTotal        = settings->GetNTotal();

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    std::cout << "PE " << settings->GetMyPE() << ": kStart is " << settings->GetKStart() << " kEnd is " << settings->GetKEnd() << std::endl;
    for( unsigned k = settings->GetKStart(); k <= settings->GetKEnd(); ++k ) {
        std::cout << "PE " << settings->GetMyPE() << ": xi = " << xi[k] << std::endl;
        uQ[k - settings->GetKStart()] = solver->Solve( xi[k] );
    }

    log->info( "" );
    log->info( "Process exited normally on PE {:03.8f} .", double( settings->GetMyPE() ) );

    // compute Moments
    MatVec u( nCells, Matrix( nStates, nTotal ) );
    MatVec uQFinalPE( nCells, Matrix( nStates, settings->GetKEnd() - settings->GetKStart() + 1 ) );
    for( unsigned j = 0; j < nCells; ++j ) {
        uQFinalPE[j].reset();
        for( unsigned l = 0; l < nStates; ++l ) {
            for( unsigned k = settings->GetKStart(); k <= settings->GetKEnd(); ++k ) {
                uQFinalPE[j]( l, k - settings->GetKStart() ) = uQ[k - settings->GetKStart()][j]( l, 0 );
            }
        }
        u[j].reset();
        multOnPENoReset( uQFinalPE[j], closure->GetPhiTildeWf(), u[j], settings->GetKStart(), settings->GetKEnd() );
    }

    auto uMoments = u;

    // perform reduction to obtain full moments on all PEs
    std::vector<int> PEforCell = settings->GetPEforCell();
    for( unsigned j = 0; j < nCells; ++j ) {
        u[j].reset();
        MPI_Reduce( uMoments[j].GetPointer(), u[j].GetPointer(), int( nStates * nTotal ), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    }

    if( settings->GetMyPE() == 0 ) {
        std::cout << u[nCells - 1] << std::endl;
        std::cout << "Exporting solution..." << std::endl;
        // export moments
        solver->Export( u );

        Matrix meanAndVar( 2 * nStates, nCells, 0.0 );
        for( unsigned j = 0; j < nCells; ++j ) {
            // expected value
            for( unsigned s = 0; s < nStates; ++s ) {
                meanAndVar( s, j ) = u[j]( s, 0 );
            }
            // variance
            for( unsigned s = 0; s < nStates; ++s ) {
                for( unsigned i = 1; i < settings->GetNTotal(); ++i ) meanAndVar( s + nStates, j ) += std::pow( u[j]( s, i ), 2 );
            }
        }

        if( settings->GetProblemType() == P_SHALLOWWATER_2D )
            mesh->ExportShallowWater( meanAndVar );
        else
            mesh->Export( meanAndVar );
    }

    delete solver;
    delete problem;
    delete mesh;
    delete settings;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
