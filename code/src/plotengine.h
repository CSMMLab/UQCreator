#ifndef PLOTENGINE_H
#define PLOTENGINE_H

#include <QWidget>
#include <blaze/math/DynamicMatrix.h>
#include <iostream>

#include "problem.h"
#include "qcustomplot.h"
#include "typedefs.h"

struct Result1D {
    Vector x_numeric;
    Vector y_numeric;
    Vector x_exact;
    Vector y_exact;
};

namespace Ui {
class PlotEngine;
}

class PlotEngine : public QWidget
{
    Q_OBJECT

  public:
    explicit PlotEngine();
    ~PlotEngine();

    void setPlot( unsigned id, Result1D data, QString title, QString xLabel, QString yLabel );

  private:
    Ui::PlotEngine* ui;
    QVector<double> BlazeToQVector( const Vector& v );

  private slots:
    void keyPressEvent( QKeyEvent* event );
};

#endif    // PLOTENGINE_H
