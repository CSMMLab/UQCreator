#include "quadrature.h"

Quadrature::Quadrature( Problem* p ) : _value( 0.0 ), _problem( p ) {
    _nodes   = _polynomial->GetNodes();
    _weights = _polynomial->GetWeights();
}

double Quadrature::Evaluate() {
    std::cerr << "[ERROR]: Not yet implemented." << std::endl;
    return -1.0;
}

Vector Quadrature::GetNodes() { return _nodes; }

Vector Quadrature::GetWeights() { return _weights; }
