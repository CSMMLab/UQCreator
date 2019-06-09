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
    sinks.push_back( terminalSink );
    if( fileLogLvl != spdlog::level::off ) sinks.push_back( fileSink );

    auto event_logger = std::make_shared<spdlog::logger>( "event", begin( sinks ), end( sinks ) );
    spdlog::register_logger( event_logger );
    spdlog::flush_every( std::chrono::seconds( 5 ) );

    auto momentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + currentDateTime() + "_moments" );
    momentFileSink->set_level( spdlog::level::info );
    momentFileSink->set_pattern( "%v" );
    auto moment_logger = std::make_shared<spdlog::logger>( "moments", momentFileSink );
    spdlog::register_logger( moment_logger );

    auto dualsFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>( outputDir + "/logs/" + currentDateTime() + "_duals" );
    dualsFileSink->set_level( spdlog::level::info );
    dualsFileSink->set_pattern( "%v" );
    auto duals_logger = std::make_shared<spdlog::logger>( "duals", dualsFileSink );
    spdlog::register_logger( duals_logger );
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

    std::cout << "PE " << settings->GetMyPE() << ": kStart is " << settings->GetKStart() << " kEnd is " << settings->GetKEnd() << std::endl;
    for( unsigned k = settings->GetKStart(); k <= settings->GetKEnd(); ++k ) {
        std::cout << "PE " << settings->GetMyPE() << ": xi = " << xi[k] << std::endl;
        uQ[k - settings->GetKStart()] = solver->Solve( xi[k] );
    }

    log->info( "" );
    log->info( "Process exited normally on PE {:03.8f} .", double( settings->GetMyPE() ) );

    // compute Moments
    MatVec u( nCells, Matrix( nStates, nTotal ) );
    MatVec uQFinal( nCells, Matrix( nStates, settings->GetKEnd() - settings->GetKStart() + 1 ) );
    for( unsigned j = 0; j < nCells; ++j ) {
        uQFinal[j].reset();
        for( unsigned l = 0; l < nStates; ++l ) {
            for( unsigned k = settings->GetKStart(); k <= settings->GetKEnd(); ++k ) {
                uQFinal[j]( l, k - settings->GetKStart() ) = uQ[k - settings->GetKStart()][j]( l, 0 );
            }
        }
        u[j].reset();
        multOnPENoReset( uQFinal[j], closure->GetPhiTildeWf(), u[j], settings->GetKStart(), settings->GetKEnd() );
    }

    auto uMoments = u;

    // perform reduction to obtain full moments on all PEs
    std::vector<int> PEforCell = settings->GetPEforCell();
    for( unsigned j = 0; j < nCells; ++j ) {
        u[j].reset();
        MPI_Reduce( uMoments[j].GetPointer(), u[j].GetPointer(), int( nStates * nTotal ), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    }

    // compute errors
    std::vector<std::vector<double>> referenceSolution;
    unsigned nQTotal = settings->GetNQTotal();
    if( settings->GetMyPE() != 0 ) return EXIT_SUCCESS;
    Matrix meanAndVar;
    Matrix meanAndVarErrors = Matrix( 2 * nStates, mesh->GetNumCells(), 0.0 );
    if( settings->HasExactSolution() ) {
        meanAndVar = Matrix( 4 * nStates, mesh->GetNumCells(), 0.0 );
    }
    else {
        meanAndVar = Matrix( 2 * nStates, mesh->GetNumCells(), 0.0 );
    }
    Matrix phiTildeWf = closure->GetPhiTildeWf();
    Vector tmp( nStates, 0.0 );
    for( unsigned j = 0; j < nCells; ++j ) {
        // expected value
        for( unsigned k = 0; k < settings->GetNQTotal(); ++k ) {

            for( unsigned i = 0; i < nStates; ++i ) {
                meanAndVar( i, j ) += uQFinal[j]( i, k ) * phiTildeWf( k, 0 );
            }
        }

        // variance
        for( unsigned k = 0; k < nQTotal; ++k ) {
            for( unsigned i = 0; i < nStates; ++i ) {
                meanAndVar( i + nStates, j ) += pow( uQFinal[j]( i, k ) - meanAndVar( i, j ), 2 ) * phiTildeWf( k, 0 );
            }
        }
    }
    // compute distance to exact solution
    if( settings->HasExactSolution() ) {
        // store xGrid on vector
        Matrix xGrid( nCells, settings->GetMeshDimension() );
        for( unsigned j = 0; j < nCells; ++j ) {
            Vector midPj = mesh->GetCenterPos( j );
            for( unsigned s = 0; s < settings->GetMeshDimension(); ++s ) {
                xGrid( j, s ) = midPj[s];
            }
        }
        Vector xiEta( settings->GetNDimXi() );
        std::vector<Polynomial*> quad( 2 );
        unsigned nQuadFine = 5 * settings->GetNQuadPoints();    // define fine quadrature for exact solution
        quad[0]            = Polynomial::Create( settings, nQuadFine, DistributionType::D_LEGENDRE );
        quad[1]            = Polynomial::Create( settings, nQuadFine, DistributionType::D_HERMITE );
        unsigned n;

        // compute indices for quad points
        std::vector<std::vector<unsigned>> indicesQ;
        unsigned nQTotal = pow( nQuadFine, settings->GetNDimXi() );
        indicesQ.resize( nQTotal );
        for( unsigned k = 0; k < nQTotal; ++k ) {
            indicesQ[k].resize( settings->GetNDimXi() );
            for( unsigned l = 0; l < settings->GetNDimXi(); ++l ) {
                indicesQ[k][l] = unsigned( ( k - k % unsigned( std::pow( nQuadFine, l ) ) ) / unsigned( std::pow( nQuadFine, l ) ) ) % nQuadFine;
            }
        }

        for( unsigned k = 0; k < nQTotal; ++k ) {
            for( unsigned l = 0; l < settings->GetNDimXi(); ++l ) {
                if( settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                xiEta[l] = quad[n]->GetNodes()[indicesQ[k][l]];
            }
            Matrix exactSolOnMesh = problem->ExactSolution( settings->GetTEnd(), xGrid, xiEta );
            for( unsigned j = 0; j < nCells; ++j ) {
                // expected value exact
                for( unsigned i = 0; i < nStates; ++i ) {
                    double wfXi = 1.0;
                    for( unsigned l = 0; l < settings->GetNDimXi(); ++l ) {
                        if( settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                        if( settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                        wfXi *= quad[n]->fXi( quad[n]->GetNodes()[indicesQ[k][l]] ) * quad[n]->GetWeights()[indicesQ[k][l]];
                    }

                    meanAndVar( 2 * nStates + i, j ) += exactSolOnMesh( j, i ) * wfXi;
                }
            }
        }
        // variance exact
        for( unsigned k = 0; k < nQTotal; ++k ) {
            for( unsigned l = 0; l < settings->GetNDimXi(); ++l ) {
                if( settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                if( settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                xiEta[l] = quad[n]->GetNodes()[indicesQ[k][l]];
            }

            Matrix exactSolOnMesh = problem->ExactSolution( settings->GetTEnd(), xGrid, xiEta );
            for( unsigned j = 0; j < nCells; ++j ) {
                // expected value exact
                for( unsigned i = 0; i < nStates; ++i ) {
                    double wfXi = 1.0;
                    for( unsigned l = 0; l < settings->GetNDimXi(); ++l ) {
                        if( settings->GetDistributionType( l ) == DistributionType::D_LEGENDRE ) n = 0;
                        if( settings->GetDistributionType( l ) == DistributionType::D_HERMITE ) n = 1;
                        wfXi *= quad[n]->fXi( quad[n]->GetNodes()[indicesQ[k][l]] ) * quad[n]->GetWeights()[indicesQ[k][l]];
                    }
                    meanAndVar( 3 * nStates + i, j ) += pow( exactSolOnMesh( j, i ) - meanAndVar( 2 * nStates + i, j ), 2 ) * wfXi;
                }
            }
        }
        // write solution on reference field
        referenceSolution.resize( nCells );
        for( unsigned j = 0; j < nCells; ++j ) {
            referenceSolution[j].resize( 2 * nStates );
            for( unsigned s = 2 * nStates; s < 4 * nStates; ++s ) {
                referenceSolution[j][s - 2 * nStates] = meanAndVar( s, j );
            }
        }
    }

    if( settings->HasReferenceFile() || settings->HasExactSolution() ) {
        Vector a( 2 );
        a[0] = -0.05;
        a[1] = -0.5;
        Vector b( 2 );
        b[0]         = 1.05;
        b[1]         = 0.5;
        auto l1Error = CalculateError( settings, mesh, referenceSolution, meanAndVar, 1, a, b );
        auto l2Error = CalculateError( settings, mesh, referenceSolution, meanAndVar, 1, a, b );
        log->info( "\nExpectation Value error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = 0; i < nStates; ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
        log->info( "\nVariance error w.r.t reference solution:" );
        log->info( "State   L1-error      L2-error" );
        for( unsigned i = nStates; i < 2 * nStates; ++i ) {
            log->info( "{:1d}       {:01.5e}   {:01.5e}", i, l1Error[i], l2Error[i] );
        }
    }

    if( settings->GetMyPE() == 0 ) {
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
