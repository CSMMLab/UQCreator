#include "HLLSolver1D.h"

HLLSolver1D::HLLSolver1D( PhysicalProblem* euler1D ) { _euler1D = euler1D; }

HLLSolver1D::~HLLSolver1D() {}

void HLLSolver1D::NumericalFlux( Vector* left, Vector* right, Vector* numericalFlux ) {
    double lamdaMin, lamdaMax;
    _euler1D->EstimateMinSpeed( left, lamdaMin );
    _euler1D->EstimateMaxSpeed( right, lamdaMax );

    if( lamdaMin >= 0 ) {
        _euler1D->AnalyticFlux( left, numericalFlux );
    }
    else if( lamdaMax <= 0 ) {
        _euler1D->AnalyticFlux( right, numericalFlux );
    }
    else {
        double factor = 1 / ( lamdaMax - lamdaMin );
        Vector fluxLeft( _euler1D->GetStateDim() ), fluxRight( _euler1D->GetStateDim() );

        _euler1D->AnalyticFlux( left, &fluxLeft );
        _euler1D->AnalyticFlux( right, &fluxRight );

        *numericalFlux = *right - *left;
        numericalFlux->Multiply( lamdaMax * lamdaMin );

        fluxLeft.Multiply( lamdaMax );
        fluxRight.Multiply( -lamdaMin );

        numericalFlux->Add( &fluxLeft );
        numericalFlux->Add( &fluxRight );

        numericalFlux->Multiply( factor );
    }
}

void HLLSolver1D::t1_NumericalFlux(
    Vector* left, Vector* t1_left, Vector* right, Vector* t1_right, Vector* numericalFlux, Vector* t1_numericalFlux ) {
    double lamdaMin, t1_lamdaMin, lamdaMax, t1_lamdaMax;
    _euler1D->t1_EstimateMinSpeed( left, t1_left, lamdaMin, t1_lamdaMin );
    _euler1D->EstimateMinSpeed( left, lamdaMin );

    _euler1D->t1_EstimateMaxSpeed( right, t1_right, lamdaMax, t1_lamdaMax );
    _euler1D->EstimateMaxSpeed( right, lamdaMax );

    if( lamdaMin >= 0 ) {
        _euler1D->t1_AnalyticFlux( left, t1_left, numericalFlux, t1_numericalFlux );
        _euler1D->AnalyticFlux( left, numericalFlux );
    }
    else if( lamdaMax <= 0 ) {
        _euler1D->t1_AnalyticFlux( right, t1_right, numericalFlux, t1_numericalFlux );
        _euler1D->AnalyticFlux( right, numericalFlux );
    }
    else {
        double t1_factor = -1 / pow( ( lamdaMax - lamdaMin ), 2 ) * t1_lamdaMax + 1 / pow( ( lamdaMax - lamdaMin ), 2 ) * t1_lamdaMin;
        double factor    = 1 / ( lamdaMax - lamdaMin );

        Vector t1_fluxLeft( _euler1D->GetStateDim() ), t1_fluxRight( _euler1D->GetStateDim() );
        Vector fluxLeft( _euler1D->GetStateDim() ), fluxRight( _euler1D->GetStateDim() );

        _euler1D->t1_AnalyticFlux( left, t1_left, &fluxLeft, &t1_fluxLeft );
        _euler1D->AnalyticFlux( left, &fluxLeft );

        _euler1D->t1_AnalyticFlux( right, t1_right, &fluxRight, &t1_fluxRight );
        _euler1D->AnalyticFlux( right, &fluxRight );

        *t1_numericalFlux = *t1_right - *t1_left;
        *numericalFlux    = *right - *left;

        t1_numericalFlux->Multiply( lamdaMax * lamdaMin );
        Vector* tmp = new Vector( _euler1D->GetStateDim() );    //!!
        *tmp        = *numericalFlux;
        tmp->Multiply( lamdaMin * t1_lamdaMax + lamdaMax * t1_lamdaMin );
        t1_numericalFlux->Add( tmp );
        numericalFlux->Multiply( lamdaMax * lamdaMin );

        t1_fluxLeft.Multiply( lamdaMax );
        *tmp = fluxLeft;
        tmp->Multiply( t1_lamdaMax );
        t1_fluxLeft.Add( tmp );
        fluxLeft.Multiply( lamdaMax );

        t1_fluxRight.Multiply( -lamdaMin );
        *tmp = fluxRight;
        tmp->Multiply( -t1_lamdaMin );
        t1_fluxRight.Add( tmp );
        fluxRight.Multiply( -lamdaMin );

        t1_numericalFlux->Add( &t1_fluxLeft );
        numericalFlux->Add( &fluxLeft );

        t1_numericalFlux->Add( &t1_fluxRight );
        numericalFlux->Add( &fluxRight );

        t1_numericalFlux->Multiply( factor );
        *tmp = *numericalFlux;
        tmp->Multiply( t1_factor );
        t1_numericalFlux->Add( tmp );
        numericalFlux->Multiply( factor );

        delete tmp;
    }
}
