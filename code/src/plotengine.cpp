#include "plotengine.h"
#include "ui_plotengine.h"

PlotEngine::PlotEngine( Settings* s, Closure* c, Mesh* m, Problem* p )
    : QWidget( nullptr ), ui( new Ui::PlotEngine ), _settings( s ), _closure( c ), _mesh( m ), _problem( p ) {
    ui->setupUi( this );
    this->setStyleSheet( "background-color:white;" );
    _plots.push_back( ui->plot0 );
    _plots.push_back( ui->plot1 );
    for( const auto& p : _plots ) {
        p->setInteractions( QCP::iRangeDrag | QCP::iRangeZoom );
    }
    this->setupPlots();

    Legendre quad( _settings->GetNQuadPoints() );
    Legendre quadFine( _nQuadFine );
    _w          = quad.GetWeights();
    _wFine      = quadFine.GetWeights();
    _xiQuad     = quad.GetNodes();
    _xiQuadFine = quadFine.GetNodes();

    _x = _mesh->GetNodePositionsX();
    _xFine.resize( _nFine );
    auto a = _mesh->GetGrid()[0]->GetNode( 0 )->coords[0];
    auto b = _mesh->GetGrid()[_mesh->GetNumCells() - 1]->GetNode( 1 )->coords[0];
    for( unsigned int j = 0; j < _nFine; ++j ) {
        _xFine[j] = a + j * ( b - a ) / ( _nFine - 1 );
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

void PlotEngine::setupPlots() {
    if( _settings->GetMeshDimension() == 1 ) {
        this->setupPlot( 0, "Mean", "x", "f(x)" );
        this->setupPlot( 1, "Var", "x", "f(x)" );
    }
    else if( _settings->GetMeshDimension() == 2 ) {
    }

    this->show();
    QCoreApplication::processEvents();
}

void PlotEngine::setupPlot( unsigned id, QString title, QString xLabel, QString yLabel ) {
    QCustomPlot* plot = _plots[id];
    plot->plotLayout()->insertRow( 0 );
    plot->plotLayout()->addElement( 0, 0, new QCPTextElement( plot, title ) );
    plot->xAxis->setLabel( xLabel );
    plot->yAxis->setLabel( yLabel );
    plot->axisRect()->setupFullAxesBox( true );
    if( _settings->GetMeshDimension() == 1 ) {
        std::vector<QColor> colors{Qt::blue, Qt::red};
        std::vector<QString> names{"numerical", "exact"};
        for( unsigned i = 0; i < 2; ++i ) {
            plot->addGraph();
            plot->graph()->setName( names[i] );
            QPen pen;
            pen.setColor( colors[i] );
            pen.setWidth( 2 );
            plot->graph()->setPen( pen );
        }
    }
    else if( _settings->GetMeshDimension() == 2 ) {
        // TODO
    }
}

void PlotEngine::updatePlotData( double time, const MatVec& lambda ) {
    ui->lTimeValue->setText( QString::fromStdString( std::to_string( time ) + " [s]" ) );
    if( _settings->GetMeshDimension() == 1 ) {
        Vector exResMean( _nFine, 0.0 ), resMean( _settings->GetNumCells(), 0.0 );

        for( unsigned int j = 0; j < _nFine; ++j ) {
            for( unsigned int k = 0; k < _nQuadFine; ++k ) {
                exResMean[j] += 0.5 * _wFine[k] * _problem->ExactSolution( _settings->GetTEnd(), _xFine[j], _xiQuadFine[k] );
            }
        }
        Vector resVec( _settings->GetNStates(), 0.0 );
        for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
            for( unsigned k = 0; k < _settings->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( lambda[j], _xiQuad, k ) );
                resMean[j] += 0.5 * _w[k] * resVec[0];
            }
        }

        Vector exResVar( _nFine, 0.0 ), resVar( _settings->GetNumCells(), 0.0 );
        for( unsigned int j = 0; j < _nFine; ++j ) {
            for( unsigned int k = 0; k < _nQuadFine; ++k ) {
                exResVar[j] += 0.5 * _wFine[k] * _problem->ExactSolution( _settings->GetTEnd(), _xFine[j], _xiQuadFine[k] );
            }
        }

        resVec.reset();
        double expectValue;
        for( unsigned j = 0; j < _settings->GetNumCells(); ++j ) {
            expectValue = 0.0;
            for( unsigned k = 0; k < _settings->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( lambda[j], _xiQuad, k ) );
                expectValue += 0.5 * _w[k] * resVec[0];
            }
            for( unsigned k = 0; k < _settings->GetNQuadPoints(); ++k ) {
                _closure->U( resVec, _closure->EvaluateLambda( lambda[j], _xiQuad, k ) );
                resVar[j] += 0.5 * _w[k] * pow( resVec[0] - expectValue, 2 );
            }
        }

        _plots[0]->graph( 0 )->setData( BlazeToQVector( _mesh->GetNodePositionsX() ), BlazeToQVector( resMean ) );
        _plots[0]->graph( 1 )->setData( BlazeToQVector( _xFine ), BlazeToQVector( exResMean ) );
        _plots[1]->graph( 0 )->setData( BlazeToQVector( _mesh->GetNodePositionsX() ), BlazeToQVector( resVar ) );
        _plots[1]->graph( 1 )->setData( BlazeToQVector( _xFine ), BlazeToQVector( exResVar ) );
    }
    this->replot();
}

/* 2D Plot
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
*/

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

void PlotEngine::refresh() { QCoreApplication::processEvents(); }
