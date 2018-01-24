#include "gnuplotlib.h"

GnuplotLib::GnuplotLib( Problem* problem ) : PlotEngine( problem ) {}

void GnuplotLib::Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 ) {
    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d( std::make_pair( x1, y1 ) ) << "with lines\n";
}

void GnuplotLib::Plot1D( const std::vector<double>& x1,
                         const std::vector<double>& y1,
                         const std::vector<double>& x2,
                         const std::vector<double>& y2 ) {
    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d( std::make_pair( x1, y1 ) ) << "with lines, " << gp.file1d( std::make_pair( x2, y2 ) ) << "with lines\n";
}

void GnuplotLib::Plot1D( const blaze::DynamicVector<double>& x, const blaze::DynamicVector<double>& y ) {
    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d( std::make_pair( this->BlazeToStdVector( x ), this->BlazeToStdVector( y ) ) ) << "with lines\n";
}

void GnuplotLib::Plot1D( const blaze::DynamicVector<double>& x1,
                         const blaze::DynamicVector<double>& y1,
                         const blaze::DynamicVector<double>& x2,
                         const blaze::DynamicVector<double>& y2 ) {
    Gnuplot gp;
    gp << "set key off\n";
    gp << "plot" << gp.file1d( std::make_pair( this->BlazeToStdVector( x1 ), this->BlazeToStdVector( y1 ) ) ) << "with lines, "
       << gp.file1d( std::make_pair( this->BlazeToStdVector( x2 ), this->BlazeToStdVector( y2 ) ) ) << "with lines\n";
}
