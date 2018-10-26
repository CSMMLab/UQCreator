#ifndef PLOTENGINE_H
#define PLOTENGINE_H

#include <QWidget>
#include <iostream>

#include "closure.h"
#include "legendre.h"
#include "mesh.h"
#include "problem.h"
#include "qcustomplot.h"
#include "settings.h"
#include "typedefs.h"

namespace Ui {
class PlotEngine;
}

class PlotEngine : public QWidget
{
    Q_OBJECT

  public:
    explicit PlotEngine( Settings* s, Closure* c, Mesh* m, Problem* p );
    ~PlotEngine();
    std::vector<QCustomPlot*> _plots;

    void setupPlots();
    void updatePlotData( double time, const MatVec& );
    void refresh();

  private:
    Ui::PlotEngine* ui;
    Settings* _settings;
    Closure* _closure;
    Mesh* _mesh;
    Problem* _problem;

    const unsigned _nQuadFine = 200;
    const unsigned _nFine     = 1000;
    Vector _w, _wFine, _xiQuad, _xiQuadFine;
    Vector _x, _xFine;

    QVector<double> BlazeToQVector( const Vector& v );
    void setupPlot( unsigned id, QString title, QString xLabel, QString yLabel );
    void replot();

  private slots:
    void keyPressEvent( QKeyEvent* event );
};

#endif    // PLOTENGINE_H
