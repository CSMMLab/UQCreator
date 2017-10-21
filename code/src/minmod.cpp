#include "minmod.h"
#include "math.h"

Minmod::Minmod(Closure* pClosure, Problem* problem): Limiter(pClosure,problem){

}

double Minmod::CalculateSlope(double u0, double u1, double u2){
    return (1/_problem->GetMesh()->GetSpacing()[0])*minmod(u2-u1,u1-u0);
}

double Minmod::minmod(double a, double b){
    double y;
    if( fabs(a) < fabs(b) && a*b > 0)
        y = a;
    else if(fabs(b) < fabs(a) && a*b >0 )
        y = b;
    else
        y = 0;
    return y;
}
