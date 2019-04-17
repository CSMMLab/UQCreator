#include <chrono>
#include <iostream>
#include <mpi.h>
#include <omp.h>

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

    auto fileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + cfilename );
    fileSink->set_level( fileLogLvl );
    fileSink->set_pattern( "%Y-%m-%d %H:%M:%S.%f | %v" );

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back( terminalSink );
    if( fileLogLvl != spdlog::level::off ) sinks.push_back( fileSink );

    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );

    auto momentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + cfilename + "_moments" );
    momentFileSink->set_level( spdlog::level::info );
    momentFileSink->set_pattern( "%v" );
    auto moment_logger = std::make_shared<spdlog::logger>( "moments", momentFileSink );
    spdlog::register_logger( moment_logger );

    auto dualsFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + cfilename + "_duals" );
    dualsFileSink->set_level( spdlog::level::info );
    dualsFileSink->set_pattern( "%v" );
    auto duals_logger = std::make_shared<spdlog::logger>( "duals", dualsFileSink );
    spdlog::register_logger( duals_logger );

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

    std::string logfilename = initLogger( spdlog::level::info, spdlog::level::info, configFile );
    MPI_Barrier( MPI_COMM_WORLD );

    auto log = spdlog::get( "event" );

    Settings* settings = new Settings( configFile );
    if( settings->GetMyPE() == 0 ) {
        PrintInit( configFile );
    }
    if( settings->HasReferenceFile() ) initErrorLogger( configFile, logfilename );
    MPI_Barrier( MPI_COMM_WORLD );
    Mesh* mesh           = Mesh::Create( settings );
    Problem* problem     = Problem::Create( settings );
    MomentSolver* solver = new MomentSolver( settings, mesh, problem );

    solver->Solve();

    log->info( "" );
    log->info( "Process exited normally on PE {:1d} .", settings->GetMyPE() );

    delete solver;
    delete problem;
    delete mesh;
    delete settings;

    MPI_Finalize();

    return EXIT_SUCCESS;
}
