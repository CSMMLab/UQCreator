#include "plotengine.h"
#include "ui_plotengine.h"

PlotEngine::PlotEngine( Problem* p ) : QWidget( nullptr ), ui( new Ui::PlotEngine ), _problem( p ) {
    ui->setupUi( this );
    this->setStyleSheet( "background-color:white;" );
    _plots.push_back( ui->plot0 );
    _plots.push_back( ui->plot1 );
    _plots.push_back( ui->plot2 );
    for( const auto& p : _plots ) {
        p->setInteractions( QCP::iRangeDrag | QCP::iRangeZoom );
    }
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

void PlotEngine::setupPlot( unsigned id, QString title, QString xLabel, QString yLabel ) {
    QCustomPlot* plot = _plots[id];
    plot->plotLayout()->insertRow( 0 );
    plot->plotLayout()->addElement( 0, 0, new QCPTextElement( plot, title ) );
    plot->xAxis->setLabel( xLabel );
    plot->yAxis->setLabel( yLabel );
    plot->axisRect()->setupFullAxesBox( true );
}

void PlotEngine::addPlotData( unsigned id, Result1D data, QString name ) {
    std::vector<QColor> colors{Qt::blue, Qt::red, Qt::green};
    QCustomPlot* plot = _plots[id];
    plot->addGraph();
    plot->graph()->addData( BlazeToQVector( data.x ), BlazeToQVector( data.y ) );
    plot->graph()->setName( name );
    QPen pen;
    pen.setColor( colors[id] );
    pen.setWidth( 2 );
    plot->graph()->setPen( pen );
}

void PlotEngine::updatePlotData( unsigned id, Result1D data, QString name ) {
    QCustomPlot* plot = _plots[id];
    int idx           = 0;
    for( int i; i < plot->graphCount(); ++i ) {
        if( plot->graph( i )->name() == name ) {
            idx = i;
        }
    }
    plot->graph( idx )->setData( BlazeToQVector( data.x ), BlazeToQVector( data.y ) );
}

void PlotEngine::addPlotData( unsigned id, Result2D data, QString name ) {
    QCustomPlot* plot = _plots[id];
    plot->axisRect()->setupFullAxesBox( true );
    QCPColorMap* colorMap = new QCPColorMap( plot->xAxis, plot->yAxis );

    int nx = static_cast<int>( data.data.columns() );
    int ny = static_cast<int>( data.data.rows() );
    colorMap->data()->setSize( nx, ny );
    colorMap->data()->setRange( QCPRange( data.xRange.first, data.xRange.second ), QCPRange( data.yRange.first, data.yRange.second ) );
    for( int xIndex = 0; xIndex < nx; ++xIndex ) {
        for( int yIndex = 0; yIndex < ny; ++yIndex ) {
            colorMap->data()->setCell( xIndex, yIndex, data.data( static_cast<unsigned>( xIndex ), static_cast<unsigned>( yIndex ) ) );
        }
    }
    QCPColorScale* colorScale = new QCPColorScale( plot );
    plot->plotLayout()->addElement( 0, 1, colorScale );
    colorScale->setType( QCPAxis::atRight );
    colorMap->setColorScale( colorScale );
    colorMap->setGradient( QCPColorGradient::gpJet );
    colorMap->setName( name );
    colorMap->rescaleDataRange();

    QCPMarginGroup* marginGroup = new QCPMarginGroup( plot );
    plot->axisRect()->setMarginGroup( QCP::msBottom | QCP::msTop, marginGroup );
    colorScale->setMarginGroup( QCP::msBottom | QCP::msTop, marginGroup );

    plot->rescaleAxes();
}

void PlotEngine::replot() {
    for( const auto& p : _plots ) {
        if( p->graphCount() == 0 ) {
            p->setVisible( false );
            continue;
        }
        p->rescaleAxes();
        p->replot();
    }
    QCoreApplication::processEvents();
}
