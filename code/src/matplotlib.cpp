#include "matplotlib.h"

Matplotlib::Matplotlib(){

}

void Matplotlib::Plot1D(const std::vector<double>& x1, const std::vector<double>& y1){
    matplotlibcpp::plot(x1,y1);
    matplotlibcpp::save("recent_plot"+_fileExt);
}

void Matplotlib::Plot1D(const std::vector<double>& x1, const std::vector<double>& y1, const std::vector<double>& x2, const std::vector<double>& y2){
    matplotlibcpp::plot(x1,y1);
    matplotlibcpp::plot(x2,y2);
    matplotlibcpp::save("recent_plot"+_fileExt);
}

void Matplotlib::Plot1D(const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1){
    matplotlibcpp::plot(this->BlazeToStdVector(x1), this->BlazeToStdVector(y1));
    matplotlibcpp::save("recent_plot"+_fileExt);
}

void Matplotlib::Plot1D(const blaze::DynamicVector<double>& x1, const blaze::DynamicVector<double>& y1, const blaze::DynamicVector<double>& x2, const blaze::DynamicVector<double>& y2){
    matplotlibcpp::plot(this->BlazeToStdVector(x1), this->BlazeToStdVector(y1));
    matplotlibcpp::plot(this->BlazeToStdVector(x2), this->BlazeToStdVector(y2));
    matplotlibcpp::save("recent_plot"+_fileExt);
}
