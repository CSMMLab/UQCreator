#include "limiter.h"
#include "minmod.h"
#include "nolimiter.h"

Limiter::Limiter( Settings* settings, Mesh* mesh, Closure* closure ) : _settings( settings ), _mesh( mesh ), _closure( closure ) {
    _dx = _mesh->GetArea( 0 );
}

Limiter::~Limiter() {}

Matrix Limiter::Slope( const Matrix& lambda1, const Matrix& lambda2, const Matrix& lambda3 ) {
    return SlopeInternal( _closure->U( _closure->EvaluateLambda( lambda1 ) ),
                          _closure->U( _closure->EvaluateLambda( lambda2 ) ),
                          _closure->U( _closure->EvaluateLambda( lambda3 ) ) );
}

double Limiter::SlopeBoundPres( const double& u, const double& slope ) {
    double M = std::max( u + 0.5 * _dx * slope, u - 0.5 * _dx * slope );
    double m = std::min( u + 0.5 * _dx * slope, u - 0.5 * _dx * slope );
    return std::min( 1.0, std::min( fabs( ( _closure->GetUPlus() - u ) / ( M - u ) ), fabs( ( _closure->GetUMinus() - u ) / ( m - u ) ) ) );
}

Matrix Limiter::SlopeInternal( const Matrix& u0, const Matrix& u1, const Matrix& u2 ) {
    unsigned nQ = _settings->GetNQuadPoints();
    double classicalSlope;
    Matrix y( _settings->GetNStates(), nQ );
    for( unsigned l = 0; l < _settings->GetNStates(); ++l ) {
        for( unsigned k = 0; k < nQ; ++k ) {
            classicalSlope = CalculateSlope( u0( l, k ), u1( l, k ), u2( l, k ) );
            y( l, k )      = classicalSlope * SlopeBoundPres( u1( l, k ), classicalSlope );
        }
    }
    return y;
}

Limiter* Limiter::Create( Settings* settings, Mesh* mesh, Closure* closure ) {
    if( settings->GetLimiterType() == LimiterType::L_MINMOD ) {
        return new Minmod( settings, mesh, closure );
    }
    else if( settings->GetLimiterType() == LimiterType::L_NONE ) {
        return new NoLimiter( settings, mesh, closure );
    }
    else {
        std::cerr << "[Limiter] Invalid limiter type!" << std::endl;
        exit( EXIT_FAILURE );
    }
}
