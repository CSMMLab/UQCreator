#include "matplotlib.h"

Matplotlib::Matplotlib( Problem* problem ) : PlotEngine( problem ) {}

Matplotlib::~Matplotlib() {}

void Matplotlib::Plot1D( const std::vector<double>& x1, const std::vector<double>& y1 ) {
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::system_clock::now().time_since_epoch() ).count();
    matplotlibcpp::plot( x1, y1 );
    matplotlibcpp::save( _outputDir + "/" + "plot_" + std::to_string( timestamp ) + _fileExt );
    matplotlibcpp::close();
}

void Matplotlib::Plot1D( const std::vector<double>& x1,
                         const std::vector<double>& y1,
                         const std::vector<double>& x2,
                         const std::vector<double>& y2 ) {
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::system_clock::now().time_since_epoch() ).count();
    matplotlibcpp::plot( x1, y1 );
    matplotlibcpp::plot( x2, y2 );
    matplotlibcpp::save( _outputDir + "/" + "plot_" + std::to_string( timestamp ) + _fileExt );
    matplotlibcpp::close();
}

void Matplotlib::Plot1D( const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1 ) {
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::system_clock::now().time_since_epoch() ).count();
    matplotlibcpp::plot( this->BlazeToStdVector( x1 ), this->BlazeToStdVector( y1 ) );
    matplotlibcpp::save( _outputDir + "/" + "plot_" + std::to_string( timestamp ) + _fileExt );
    matplotlibcpp::close();
}

void Matplotlib::Plot1D( const blaze::DynamicVector<double>& x1,
                         const blaze::DynamicVector<double>& y1,
                         const blaze::DynamicVector<double>& x2,
                         const blaze::DynamicVector<double>& y2 ) {
    auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::system_clock::now().time_since_epoch() ).count();
    matplotlibcpp::plot( this->BlazeToStdVector( x1 ), this->BlazeToStdVector( y1 ) );
    matplotlibcpp::plot( this->BlazeToStdVector( x2 ), this->BlazeToStdVector( y2 ) );
    matplotlibcpp::save( _outputDir + "/" + "plot_" + std::to_string( timestamp ) + _fileExt );
    matplotlibcpp::close();
}
