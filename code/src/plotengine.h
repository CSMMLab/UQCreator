#ifndef PLOTENGINE_H
#define PLOTENGINE_H

#include <QWidget>
#include <blaze/math/DynamicMatrix.h>
#include <iostream>

#include "problem.h"
#include "qcustomplot.h"
#include "typedefs.h"

struct Result1D {
    Vector x;
    Vector y;
};

struct Result2D {
    std::pair<double, double> xRange;
    std::pair<double, double> yRange;
    Matrix data;
};

namespace Ui {
class PlotEngine;
}

class PlotEngine : public QWidget
{
    Q_OBJECT

  public:
    explicit PlotEngine( Problem* p );
    ~PlotEngine();
    std::vector<QCustomPlot*> _plots;

    void setupPlot( unsigned id, QString title, QString xLabel, QString yLabel );
    void addPlotData( unsigned id, Result1D data, QString name );
    void addPlotData( unsigned id, Result2D data, QString name );
    void updatePlotData( unsigned id, Result1D data, QString name );
    void replot();

  private:
    Ui::PlotEngine* ui;
    Problem* _problem;
    QVector<double> BlazeToQVector( const Vector& v );

  private slots:
    void keyPressEvent( QKeyEvent* event );
};

#endif    // PLOTENGINE_H
