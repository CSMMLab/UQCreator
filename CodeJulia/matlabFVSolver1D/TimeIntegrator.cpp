#include "TimeIntegrator.h"

TimeIntegrator::TimeIntegrator( Settings* settings, Mesh* mesh, Rhs* rhs ) : _settings( settings ), _mesh( mesh ), _rhs( rhs ) {}
TimeIntegrator::~TimeIntegrator() {}
