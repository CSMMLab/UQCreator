#include "limiter.h"

Limiter::Limiter(Closure* pClosure):_closure(pClosure){

}

blaze::DynamicVector<double> Limiter::Slope(const blaze::DynamicVector<double>& lambda1, const blaze::DynamicVector<double>& lambda2, const blaze::DynamicVector<double>& lambda3){
    return SlopeInternal(_closure->UKinetic(_closure->EvaluateLambda(lambda1)),_closure->UKinetic(_closure->EvaluateLambda(lambda2)),_closure->UKinetic(_closure->EvaluateLambda(lambda3)));
}

blaze::DynamicVector<double> SlopeInternal(blaze::DynamicVector<double> u0, blaze::DynamicVector<double> u1, blaze::DynamicVector<double> u2){
    int nQ = u0.size();
    blaze::DynamicVector<double> y(nQ);
    for( int k = 0; k<nQ; ++k ){
        y[k] = CalculateSlope(u0[k],u1[k],u2[k])*slopeBoundPreserving(u0[k],u1[k],u2[k]);
    }
}
