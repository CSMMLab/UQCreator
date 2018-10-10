#include "Mesh1D.h"

Mesh1D::Mesh1D( double xLeft, double xRight, unsigned int dimX, PhysicalProblem* physProb, Settings* settings )
    : Mesh( dimX, settings ), _xLeft( xLeft ), _xRight( xRight ), _dimX( dimX ) {
    // set global variables
    _lengthX = _xRight - _xLeft;
    _dx      = _lengthX / ( _dimX );
    _cells.reserve( _dimX );

    // fill mesh with 1D Cells
    for( unsigned int i = 0; i < _dimX; i++ ) {
        Cell* cell = new Cell1D( physProb->GetStateDim() );
        _cells[i]  = cell;
    }
}

Mesh1D::~Mesh1D() {
    for( unsigned int i = 0; i < _cells.size(); i++ ) delete _cells[i];
}

Cell* Mesh1D::GetCell( int i ) const {
    if( i >= 0 && i < _nCells ) {
        return _cells[i];
    }
    else if( i < 0 && _settings->GetBounCon() != PERIODIC ) {
        return _cells[0];
    }
    else if( i > _nCells - 1 && _settings->GetBounCon() != PERIODIC ) {
        return _cells[_nCells - 1];
    }
    else if( i < 0 && _settings->GetBounCon() == PERIODIC ) {
        if( i == -1 ) return _cells[_nCells - 1];
        if( i == -2 ) return _cells[_nCells - 2];
    }
    else if( i > _nCells - 1 && _settings->GetBounCon() == PERIODIC ) {
        if( i == _nCells ) return _cells[0];
        if( i == _nCells + 1 ) return _cells[1];
    }
    return NULL;
}

Cell* Mesh1D::GetLeftCell( int i ) const {
    if( i > _nCells - 1 && _settings->GetBounCon() != PERIODIC ) {
        return _cells[_nCells - 1];
    }
    if( i <= 0 && _settings->GetBounCon() == OPEN ) return _cells[0];
    if( i <= 0 && _settings->GetBounCon() == PERIODIC ) {
        // std::cout<<"GetLeftCell called by cell "<<i<<std::endl;
        if( _settings->GetLimiter() == NOLIMITER ) {
            return _cells[_nCells - 1];
        }
        if( _settings->GetLimiter() == MINMOD || _settings->GetLimiter() == VANLEER ) {
            if( i == 0 )
                return _cells[_nCells - 1];
            else
                return _cells[_nCells - 2];
        }
    }
    if( i > _nCells - 1 && _settings->GetBounCon() == PERIODIC ) {
        if( _settings->GetLimiter() == NOLIMITER ) {
            return _cells[_nCells - 1];
        }
        if( _settings->GetLimiter() == MINMOD || _settings->GetLimiter() == VANLEER ) {
            if( i == _nCells )
                return _cells[_nCells - 1];
            else
                return _cells[0];
        }
    }
    // if( i <= 0 && _settings->GetBounCon() == DIRICHLET){
    // ghostCellLeft u,rho = cell[0], p = dirichlet
    // return _ghostCellLeft;
    //}
    return _cells[i - 1];
}

Cell* Mesh1D::GetRightCell( int i ) const {
    if( i < 0 && _settings->GetBounCon() != PERIODIC ) {
        return _cells[0];
    }
    if( i >= _nCells - 1 && _settings->GetBounCon() == OPEN ) return _cells[_nCells - 1];
    if( i >= _nCells - 1 && _settings->GetBounCon() == PERIODIC ) {
        if( _settings->GetLimiter() == NOLIMITER ) {
            return _cells[0];
        }
        if( _settings->GetLimiter() == MINMOD || _settings->GetLimiter() == VANLEER ) {
            if( i == _nCells - 1 )
                return _cells[0];
            else
                return _cells[1];
        }
    }
    if( i < 0 && _settings->GetBounCon() == PERIODIC ) {
        if( _settings->GetLimiter() == NOLIMITER ) {
            return _cells[_nCells - 1];
        }
        if( _settings->GetLimiter() == MINMOD || _settings->GetLimiter() == VANLEER ) {
            if( i == -1 )
                return _cells[0];
            else
                return _cells[_nCells - 1];
        }
    }
    // if( i <= 0 && _settings->GetBounCon() == DIRICHLET){
    // ghostCellLeft u,rho = cell[0], p = dirichlet
    // return _ghostCellLeft;
    //}
    return _cells[i + 1];
}

double Mesh1D::GetDx() const { return _dx; }

double Mesh1D::GetXLeft() const { return _xLeft; }

double Mesh1D::GetXRight() const { return _xRight; }

int Mesh1D::GetDimX() const { return _dimX; }
