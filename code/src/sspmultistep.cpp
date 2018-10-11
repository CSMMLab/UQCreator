#include "sspmultistep.h"

SSPMultiStep::SSPMultiStep( Problem* problem, Closure* closure ) : TimeSolver( problem ) {
    _counter = 0;
    _heun    = new Heun( problem, closure );
}

SSPMultiStep::~SSPMultiStep() { delete _heun; }
