#include "plotengine.h"
#include "gnuplotlib.h"
#include "matplotlib.h"
#include "noplot.h"

PlotEngine::PlotEngine( Problem* problem ) : _problem( problem ) { _outputDir = _problem->GetOutputDir(); }

PlotEngine::~PlotEngine() {}

PlotEngine* PlotEngine::Create( Problem* problem ) {
    auto file          = cpptoml::parse_file( problem->GetInputFile() );
    auto general       = file->get_table( "plot" );
    std::string engine = general->get_as<std::string>( "engine" ).value_or( "" );
    if( engine.compare( "gnuplot" ) == 0 ) {
        return new GnuplotLib( problem );
    }
    else if( engine.compare( "matplotlib" ) == 0 ) {
        return new Matplotlib( problem );
    }
    else if( engine.compare( "off" ) == 0 ) {
        return new noPlot( problem );
    }
    else {
        std::cerr << "Invalid plot engine type" << std::endl;
        exit( EXIT_FAILURE );
        return NULL;
    }
}

std::vector<double> PlotEngine::BlazeToStdVector( const blaze::DynamicVector<double>& v ) {
    std::vector<double> ret( v.size(), 0.0 );
    for( unsigned i = 0; i < v.size(); ++i ) {
        ret[i] = v[i];
    }
    return ret;
}
