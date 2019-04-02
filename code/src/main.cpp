#include <chrono>
#include <iostream>
#include <mpi.h>
#include <omp.h>

#include "closure.h"
#include "mesh.h"
#include "momentsolver.h"
#include "problem.h"
#include "settings.h"

#include "spdlog/spdlog.h"

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_sinks.h"

bool CheckInput( std::string& configFile, int argc, char* argv[] ) {
    std::string usage_help = "\n"
                             "Usage: " +
                             std::string( argv[0] ) + " -c inputfile\n";

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
    std::string filename = currentDateTime();
    int ctr              = 0;
    if( std::filesystem::exists( outputDir + "/logs/" + filename ) ) {
        filename += "_" + std::to_string( ++ctr );
    }
    while( std::filesystem::exists( outputDir + "/logs/" + filename ) ) {
        filename.pop_back();
        filename += std::to_string( ++ctr );
    }

    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename );
    fileSink->set_level( fileLogLvl );
    fileSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back( terminalSink );
    if( fileLogLvl != spdlog::level::off ) sinks.push_back( fileSink );

    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );

    auto momentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_moments" );
    momentFileSink->set_level( spdlog::level::info );
    momentFileSink->set_pattern( "%v" );
    auto moment_logger = std::make_shared<spdlog::logger>( "moments", momentFileSink );
    spdlog::register_logger( moment_logger );

    auto dualsFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + filename + "_duals" );
    dualsFileSink->set_level( spdlog::level::info );
    dualsFileSink->set_pattern( "%v" );
    auto duals_logger = std::make_shared<spdlog::logger>( "duals", dualsFileSink );
    spdlog::register_logger( duals_logger );
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
            if( line[0] != '#' ) log->info( " {0}", line );
        }
    }
    log->info( "================================================================" );
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

    Settings* settings = new Settings( configFile );
    if( settings->GetMyPE() == 0 ) PrintInit( configFile );

    // DEBUG START
    Problem* p = Problem::Create( settings );
    Closure* c = Closure::Create( settings );

    Matrix u( settings->GetNStates(), settings->GetNTotal(), 0.0 );
    Matrix lambda( settings->GetNStates(), settings->GetNTotal(), 0.0 );

    u( 0, 0 ) = 0.95452;
    u( 0, 1 ) = -0.17321;
    u( 0, 2 ) = 0.0493441;
    u( 0, 3 ) = 0.054137;
    u( 0, 4 ) = -0.070956;

    u( 1, 0 ) = 324.193;
    u( 1, 1 ) = -3.92017;
    u( 1, 2 ) = -0.378675;
    u( 1, 3 ) = 0.470835;
    u( 1, 4 ) = 1.47966;

    u( 2, 0 ) = -26.3248;
    u( 2, 1 ) = -0.232118;
    u( 2, 2 ) = -0.042907;
    u( 2, 3 ) = -0.254564;
    u( 2, 4 ) = 1.42992;

    u( 3, 0 ) = 228267.0;
    u( 3, 1 ) = -36289.7;
    u( 3, 2 ) = 10485.8;
    u( 3, 3 ) = 11206.6;
    u( 3, 4 ) = -15168.4;

    Vector ds( settings->GetNStates() );
    Vector u0( settings->GetNStates() );
    for( unsigned l = 0; l < settings->GetNStates(); ++l ) {
        u0[l] = u( l, 0 );
    }
    c->DS( ds, u0 );
    for( unsigned l = 0; l < settings->GetNStates(); ++l ) {
        lambda( l, 0 ) = ds[l];
    }
    std::cout << "initial lambda " << lambda << std::endl;
    c->SolveClosure( lambda, u, settings->GetNTotal() - 1, settings->GetNQTotal() );
    std::cout << "closure solved with " << lambda << std::endl;

    MPI_Finalize();

    return EXIT_SUCCESS;
    // DEBUG END

    // Mesh* mesh           = Mesh::Create( settings );
    // Problem* problem     = Problem::Create( settings );
    // MomentSolver* solver = new MomentSolver( settings, mesh, problem );

    // solver->Solve();

    log->info( "" );
    log->info( "Process exited normally on PE {:1d} .", settings->GetMyPE() );

    // delete solver;
    // delete problem;
    // delete mesh;
    delete settings;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
