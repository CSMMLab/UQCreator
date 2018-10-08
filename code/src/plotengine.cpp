#include "plotengine.h"
#include "ui_plotengine.h"

PlotEngine::PlotEngine() : QWidget( nullptr ), ui( new Ui::PlotEngine ) {
    ui->setupUi( this );
    this->setStyleSheet( "background-color:white;" );
    ui->plot0->setInteractions( QCP::iRangeDrag | QCP::iRangeZoom );
    ui->plot1->setInteractions( QCP::iRangeDrag | QCP::iRangeZoom );
    ui->plot2->setInteractions( QCP::iRangeDrag | QCP::iRangeZoom );
}

PlotEngine::~PlotEngine() { delete ui; }

QVector<double> PlotEngine::BlazeToQVector( const Vector& v ) {
    QVector<double> ret( static_cast<int>( v.size() ), 0.0 );
    for( unsigned i = 0; i < v.size(); ++i ) {
        ret[static_cast<int>( i )] = v[i];
    }
    return ret;
}

void PlotEngine::keyPressEvent( QKeyEvent* event ) {
    auto keyCode = event->key();
    if( keyCode == 83 ) {
        /*
        QString filename = QFileDialog::getSaveFileName( nullptr, tr( "Save plot as" ), QDir::homePath(), tr( "Vector graphic (*.pdf)" ) );
        if( filename.compare( "" ) != 0 ) {
            if( QFileInfo( filename ).suffix().compare( "pdf" ) != 0 ) {
                filename += ".pdf";
            }
            ui->plot0->savePdf( filename );
        }
        */
    }
}

void PlotEngine::setPlot( unsigned id, Result1D data, QString title, QString xLabel, QString yLabel ) {
    QCustomPlot* plot = nullptr;
    if( id == 0 ) {
        plot = ui->plot0;
    }
    else if( id == 1 ) {
        plot = ui->plot1;
    }
    else if( id == 2 ) {
        plot = ui->plot2;
    }

    plot->addGraph();
    plot->graph()->addData( BlazeToQVector( data.x_numeric ), BlazeToQVector( data.y_numeric ) );
    QPen pen;
    pen.setColor( Qt::blue );
    pen.setWidth( 2 );
    plot->graph()->setPen( pen );
    plot->addGraph();
    plot->graph()->addData( BlazeToQVector( data.x_exact ), BlazeToQVector( data.y_exact ) );
    pen.setColor( Qt::red );
    pen.setStyle( Qt::DashLine );
    plot->graph()->setPen( pen );
    plot->axisRect()->setupFullAxesBox( true );
    plot->xAxis->setLabel( xLabel );
    plot->yAxis->setLabel( yLabel );
    plot->plotLayout()->insertRow( 0 );
    plot->plotLayout()->addElement( 0, 0, new QCPTextElement( plot, title ) );
    plot->rescaleAxes();
}
