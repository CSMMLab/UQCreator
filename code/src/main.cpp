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

Vector CalculateError( const Matrix& solution, Settings* settings, Mesh* mesh, unsigned LNorm, const Vector& a, const Vector& b ) {
    auto referenceSolution = mesh->Import();
    Vector error( 2 * settings->GetNStates(), 0.0 );
    Vector refNorm( 2 * settings->GetNStates(), 0.0 );
    if( settings->GetMeshDimension() == 2 ) {
        for( unsigned j = 0; j < settings->GetNumCells(); ++j ) {
            if( mesh->GetGrid()[j]->GetCenter()[0] > a[0] && mesh->GetGrid()[j]->GetCenter()[0] < b[0] && mesh->GetGrid()[j]->GetCenter()[1] > a[1] &&
                mesh->GetGrid()[j]->GetCenter()[1] < b[1] ) {

                switch( LNorm ) {
                    case 0:
                        for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
                            error[s]   = std::max( std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j ), error[s] );
                            refNorm[s] = std::max( std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j ), refNorm[s] );
                        }
                        break;
                    case 1:
                        for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
                            error[s] += std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j );
                            refNorm[s] += std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j );
                        }
                        break;
                    case 2:
                        for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
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
        for( unsigned j = 0; j < settings->GetNumCells(); ++j ) {
            switch( LNorm ) {
                case 0:
                    for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
                        error[s]   = std::max( std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j ), error[s] );
                        refNorm[s] = std::max( std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j ), refNorm[s] );
                    }
                    break;
                case 1:
                    for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
                        error[s] += std::fabs( ( solution( s, j ) - referenceSolution[j][s] ) ) * mesh->GetArea( j );
                        refNorm[s] += std::fabs( referenceSolution[j][s] ) * mesh->GetArea( j );
                    }
                    break;
                case 2:
                    for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
                        error[s] += std::pow( ( solution( s, j ) - referenceSolution[j][s] ), 2 ) * mesh->GetArea( j );
                        refNorm[s] += std::pow( referenceSolution[j][s], 2 ) * mesh->GetArea( j );
                    }

                    break;
                default: exit( EXIT_FAILURE );
            }
        }
    }
    for( unsigned s = 0; s < 2 * settings->GetNStates(); ++s ) {
        error[s] = std::pow( error[s] / refNorm[s], 1.0 / std::max( double( LNorm ), 1.0 ) );
    }
    return error;
}

void WriteErrors( const MatVec& u, Settings* settings, Mesh* mesh ) {
    Matrix meanAndVar( 2 * settings->GetNStates(), settings->GetNumCells(), 0.0 );
    for( unsigned j = 0; j < settings->GetNumCells(); ++j ) {
        // expected value
        for( unsigned s = 0; s < settings->GetNStates(); ++s ) {
            meanAndVar( s, j ) = u[j]( s, 0 );
        }
        // variance
        for( unsigned s = 0; s < settings->GetNStates(); ++s ) {
            for( unsigned i = 1; i < settings->GetNTotal(); ++i ) meanAndVar( s + settings->GetNStates(), j ) += std::pow( u[j]( s, i ), 2 );
        }
    }

    if( settings->GetProblemType() == P_SHALLOWWATER_2D )
        mesh->ExportShallowWater( meanAndVar );
    else
        mesh->Export( meanAndVar );

    if( settings->HasReferenceFile() ) {
        // define rectangle for error computation
        Vector a( 2 );
        a[0] = -0.05;
        a[1] = -0.5;
        Vector b( 2 );
        b[0]                  = 1.05;
        b[1]                  = 0.5;
        auto l1ErrorMeanLog   = spdlog::get( "l1ErrorMean" );
        auto l2ErrorMeanLog   = spdlog::get( "l2ErrorMean" );
        auto lInfErrorMeanLog = spdlog::get( "lInfErrorMean" );
        auto l1ErrorVarLog    = spdlog::get( "l1ErrorVar" );
        auto l2ErrorVarLog    = spdlog::get( "l2ErrorVar" );
        auto lInfErrorVarLog  = spdlog::get( "lInfErrorVar" );

        auto l1Error   = CalculateError( meanAndVar, settings, mesh, 1, a, b );
        auto l2Error   = CalculateError( meanAndVar, settings, mesh, 2, a, b );
        auto lInfError = CalculateError( meanAndVar, settings, mesh, 0, a, b );

        std::ostringstream osL1ErrorMean, osL2ErrorMean, osLInfErrorMean, osL1ErrorVar, osL2ErrorVar, osLInfErrorVar;
        for( unsigned i = 0; i < settings->GetNStates(); ++i ) {
            osL1ErrorMean << std::scientific << l1Error[i] << "\t";
            osL2ErrorMean << std::scientific << l2Error[i] << "\t";
            osLInfErrorMean << std::scientific << lInfError[i] << "\t";
            osL1ErrorVar << std::scientific << l1Error[i + settings->GetNStates()] << "\t";
            osL2ErrorVar << std::scientific << l2Error[i + settings->GetNStates()] << "\t";
            osLInfErrorVar << std::scientific << lInfError[i + settings->GetNStates()] << "\t";
        }

        l1ErrorMeanLog->info( osL1ErrorMean.str() );
        l2ErrorMeanLog->info( osL2ErrorMean.str() );
        lInfErrorMeanLog->info( osLInfErrorMean.str() );
        l1ErrorVarLog->info( osL1ErrorVar.str() );
        l2ErrorVarLog->info( osL2ErrorVar.str() );
        lInfErrorVarLog->info( osLInfErrorVar.str() );

        auto log = spdlog::get( "event" );
        log->info( "\nExpectation Value error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = 0; i < settings->GetNStates(); ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
        log->info( "\nVariance error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = settings->GetNStates(); i < 2 * settings->GetNStates(); ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
    }
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

std::string initLogger( spdlog::level::level_enum terminalLogLvl, spdlog::level::level_enum fileLogLvl, std::string configFile ) {
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

    int pe;
    MPI_Comm_rank( MPI_COMM_WORLD, &pe );
    char cfilename[1024];
    if( pe == 0 ) {
        std::string filename = currentDateTime();
        int ctr              = 0;
        if( std::filesystem::exists( outputDir + "/logs/" + filename ) ) {
            filename += "_" + std::to_string( ++ctr );
        }
        while( std::filesystem::exists( outputDir + "/logs/" + filename ) ) {
            filename.pop_back();
            filename += std::to_string( ++ctr );
        }
        strncpy( cfilename, filename.c_str(), sizeof( cfilename ) );
        cfilename[sizeof( cfilename ) - 1] = 0;
    }
    MPI_Bcast( &cfilename, sizeof( cfilename ), MPI_CHAR, 0, MPI_COMM_WORLD );
    MPI_Barrier( MPI_COMM_WORLD );

    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + currentDateTime() + "_PE" + std::to_string( pe ) );
    fileSink->set_level( fileLogLvl );
    fileSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );

    std::vector<spdlog::sink_ptr> sinks;
    if( pe == 0 ) sinks.push_back( terminalSink );
    if( fileLogLvl != spdlog::level::off ) sinks.push_back( fileSink );

    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );

    auto momentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + cfilename + "_moments" );
    momentFileSink->set_level( spdlog::level::info );
    momentFileSink->set_pattern( "%v" );
    auto moment_logger = std::make_shared<spdlog::logger>( "moments", momentFileSink );
    spdlog::register_logger( moment_logger );

    return cfilename;
}

void initErrorLogger( std::string configFile, std::string filename ) {
    auto file      = cpptoml::parse_file( configFile );
    auto general   = file->get_table( "general" );
    auto outputDir = general->get_as<std::string>( "outputDir" ).value_or( "." );

    auto l1ErrorMeanSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_L1ErrorMean" );
    l1ErrorMeanSink->set_level( spdlog::level::info );
    l1ErrorMeanSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto l1ErrorMeanLogger = std::make_shared<spdlog::logger>( "l1ErrorMean", l1ErrorMeanSink );
    spdlog::register_logger( l1ErrorMeanLogger );

    auto l2ErrorMeanSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_L2ErrorMean" );
    l2ErrorMeanSink->set_level( spdlog::level::info );
    l2ErrorMeanSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto l2ErrorMeanLogger = std::make_shared<spdlog::logger>( "l2ErrorMean", l2ErrorMeanSink );
    spdlog::register_logger( l2ErrorMeanLogger );

    auto lInfErrorMeanSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_LInfErrorMean" );
    lInfErrorMeanSink->set_level( spdlog::level::info );
    lInfErrorMeanSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto lInfErrorMeanLogger = std::make_shared<spdlog::logger>( "lInfErrorMean", lInfErrorMeanSink );
    spdlog::register_logger( lInfErrorMeanLogger );

    auto l1ErrorVarSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_L1ErrorVar" );
    l1ErrorVarSink->set_level( spdlog::level::info );
    l1ErrorVarSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto l1ErrorVarLogger = std::make_shared<spdlog::logger>( "l1ErrorVar", l1ErrorVarSink );
    spdlog::register_logger( l1ErrorVarLogger );

    auto l2ErrorVarSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_L2ErrorVar" );
    l2ErrorVarSink->set_level( spdlog::level::info );
    l2ErrorVarSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto l2ErrorVarLogger = std::make_shared<spdlog::logger>( "l2ErrorVar", l2ErrorVarSink );
    spdlog::register_logger( l2ErrorVarLogger );

    auto lInfErrorVarSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_LInfErrorVar" );
    lInfErrorVarSink->set_level( spdlog::level::info );
    lInfErrorVarSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );
    auto lInfErrorVarLogger = std::make_shared<spdlog::logger>( "lInfErrorVar", lInfErrorVarSink );
    spdlog::register_logger( lInfErrorVarLogger );
}

void PrintInit( std::string configFile ) {
    auto log = spdlog::get( "event" );
    log->info( "UQCreator" );
    log->info( "================================================================" );
    log->info( "Git commit hash:\t{0}", GIT_HASH );
    log->info( "Config file:\t{0}", configFile );
    int nprocs;
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    log->info( "MPI Threads:\t{0}", nprocs );
    log->info( "OMP Threads:\t{0}", omp_get_max_threads() );
    log->info( "================================================================" );
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

    // initLogger( spdlog::level::info, spdlog::level::info, configFile );
    std::string logfilename = initLogger( spdlog::level::info, spdlog::level::info, configFile );

    auto log = spdlog::get( "event" );

    // PrintInit( configFile );

    Settings* settings = new Settings( configFile );
    if( settings->HasReferenceFile() ) initErrorLogger( configFile, logfilename );
    QuadratureGrid* quad = QuadratureGrid::Create( settings );
    settings->SetNQTotal( quad->GetNodeCount() );
    if( settings->GetMyPE() == 0 ) PrintInit( configFile );
    Mesh* mesh           = Mesh::Create( settings );
    Problem* problem     = Problem::Create( settings );
    MomentSolver* solver = new MomentSolver( settings, mesh, problem );
    Closure* closure     = Closure::Create( settings, quad );

    std::vector<MatVec> uQ;
    uQ.resize( settings->GetKEnd() - settings->GetKStart() + 1 );
    std::vector<Vector> xi = quad->GetNodes();
    unsigned nCells        = settings->GetNumCells();
    unsigned nStates       = settings->GetNStates();
    unsigned nTotal        = settings->GetNTotal();

    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();

    // std::cout << "PE " << settings->GetMyPE() << ": kStart is " << settings->GetKStart() << " kEnd is " << settings->GetKEnd() << std::endl;
    for( unsigned k = settings->GetKStart(); k <= settings->GetKEnd(); ++k ) {
        // std::cout << "PE " << settings->GetMyPE() << ": xi = " << xi[k] << std::endl;
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
        std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
        log->info( "" );
        log->info( "Finished!" );
        log->info( "" );
        log->info( "Runtime: {0}s", std::chrono::duration_cast<std::chrono::milliseconds>( toc - tic ).count() / 1000.0 );

        // std::cout << u[nCells - 1] << std::endl;
        // std::cout << "Exporting solution..." << std::endl;
        // export moments
        solver->Export( u );
        WriteErrors( u, settings, mesh );
    }

    delete solver;
    delete problem;
    delete mesh;
    delete settings;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
