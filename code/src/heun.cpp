#include "heun.h"

Heun::Heun( Problem* _problem, Closure* closure ) : TimeSolver( _problem ), _closure( closure ) {}

Heun::~Heun() {}
