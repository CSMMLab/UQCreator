#include "limiter.h"
#include <algorithm>

Limiter::Limiter(Closure* pClosure, Problem* problem):_closure(pClosure), _problem(problem) {

}

blaze::DynamicVector<double> Limiter::Slope(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2, const blaze::DynamicVector<double>& lambda3){
    return SlopeInternal(_closure->UKinetic(_closure->EvaluateLambda(lambda1)),_closure->UKinetic(_closure->EvaluateLambda(lambda2)),_closure->UKinetic(_closure->EvaluateLambda(lambda3)));
}

double Limiter::SlopeBoundPres(double u, double slope){
    double dx = _problem->GetMesh()->GetSpacing()[0];
    double M = std::max(u+0.5*dx*slope,u-0.5*dx*slope);
    double m = std::min(u+0.5*dx*slope,u-0.5*dx*slope);
    return std::min( 1.0, std::min(fabs( (_closure->GetUPlus()-u)/(M-u) ), fabs( (_closure->GetUMinus()-u)/(m-u))) );
}

blaze::DynamicVector<double> Limiter::SlopeInternal(blaze::DynamicVector<double> u0, blaze::DynamicVector<double> u1, blaze::DynamicVector<double> u2){
    int nQ = u0.size();
    double classicalSlope;
    blaze::DynamicVector<double> y(nQ);
    for( int k = 0; k<nQ; ++k ){
        classicalSlope = CalculateSlope(u0[k],u1[k],u2[k]);
        y[k] = classicalSlope*SlopeBoundPres(u1[k],classicalSlope);
    }
    return y;
}
