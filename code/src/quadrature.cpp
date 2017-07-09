#include "quadrature.h"

Quadrature::Quadrature(Problem* p) : _value(0.0), _problem(p){
    /*
    if( _problem->GetQuadType() == QUAD_TYPE_LEGENDRE ){
        _polynomial = new Legendre(_problem->GetNQuadPoints());
    }else if( _problem->GetQuadType() == QUAD_TYPE_HERMITE ){
        _polynomial = new Hermite(_problem->GetNQuadPoints());
    }else{
        std::cerr<<"[ERROR]: Quadrature unknown."<<std::endl;
    }*/
    _nodes = _polynomial->getNodes();
    _weights = _polynomial->getWeights();
}

double Quadrature::evaluate(){
    std::cerr<<"[ERROR]: Not yet implemented."<<std::endl;
    return -1.0;
}

vector Quadrature::getNodes(){
    return _nodes;
}

vector Quadrature::getWeights(){
    return _weights;
}
