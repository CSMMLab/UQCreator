#include "Cell.h"

Cell::Cell( int dim )    //: _north(NULL), _east(NULL), _south(NULL), _west(NULL)
{
    _state = new Vector( dim );
}
Cell::~Cell() { delete _state; }
Vector* Cell::GetState() { return _state; }
void Cell::SetState( Vector* vec ) { *( _state ) = *( vec ); }
/*
Cell* North()const{
    return _north;
}

Cell* East()const{
    return _east;
}

Cell* South()const{
    return _south;
}

Cell* West()const{
    return _west;
}

void SetNorth(Cell* in){
    _north = in;
}

void SetEast(Cell* in){
    _east = in;
}

void SetSouth(Cell* in){
    _south = in;
}

void SetWest(Cell* in){
    _west = in;
}*/
