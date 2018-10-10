#include "Rhs.h"

Rhs::Rhs( Settings* settings, PhysicalProblem* physProb, double dx ) : _settings( settings ), _physProb( physProb ), _dx( dx ) {}

Rhs::~Rhs() {
    // _numericalFlux is deleated by Child classes
}
