#include "plotengine.h"
#include "gnuplotlib.h"
//#include "matplotlib.h"

PlotEngine::PlotEngine() {}

PlotEngine::~PlotEngine() {}

PlotEngine *PlotEngine::Create( std::string inputFile ) {
    auto file = cpptoml::parse_file( inputFile );
    auto general = file->get_table( "plot" );
    std::string timestepping = general->get_as<std::string>( "engine" ).value_or( "" );
    if( timestepping.compare( "gnuplot" ) == 0 ) {
        return new GnuplotLib();
    }
    /*
        else if( timestepping.compare( "matplotlib" ) == 0 ) {
            return new Matplotlib();
        }
    */
    else {
        std::cerr << "Invalid plot engine type" << std::endl;
        exit( EXIT_FAILURE );
        return NULL;
    }
}

std::vector<double> PlotEngine::BlazeToStdVector( const blaze::DynamicVector<double> &v ) {
    std::vector<double> ret( v.size(), 0.0 );
    for( unsigned i = 0; i < v.size(); ++i ) {
        ret[i] = v[i];
    }
    return ret;
}
