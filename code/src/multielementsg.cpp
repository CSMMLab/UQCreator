#include "multielementsg.h"

MultiElementSG::MultiElementSG( Settings* settings ) : Closure( settings ) { _alpha = 1.0; }

MultiElementSG::~MultiElementSG() {}

void MultiElementSG::U( Vector& out, const Vector& Lambda ) { out = Lambda; }

void MultiElementSG::U( Matrix& out, const Matrix& Lambda ) { out = Lambda; }

Matrix MultiElementSG::U( const Matrix& Lambda ) { return Lambda; }

void MultiElementSG::DU( Matrix& y, const Vector& Lambda ) { y = VectorSpace::IdentityMatrix<double>( _nStates ); }

void MultiElementSG::SolveClosure( Matrix& lambda, const Matrix& u, unsigned refLevel ) { lambda = u; }

void MultiElementSG::SolveClosureSafe( Matrix& lambda, const Matrix& u, unsigned refLevel ) { lambda = u; }

Vector MultiElementSG::EvaluateLambda( const Matrix& lambda, unsigned k, unsigned nTotal ) {
    Vector out( _nStates, 0.0 );
    for( unsigned s = 0; s < _nStates; ++s ) {
        for( unsigned i = 0; i < nTotal; ++i ) {
            out[s] += lambda( s, i ) * _phiTildeVec[k][i];
        }
    }
    return out;
}

Matrix MultiElementSG::EvaluateLambda( const Matrix& lambda ) const {
    // for( unsigned n = 0; n < _settings->GetNMultiElements(); ++n ) {
    //    for( k = _set ) }

    return lambda * _phiTildeTrans;
}

Matrix MultiElementSG::EvaluateLambdaOnPE( const Matrix& lambda, unsigned levelOld, unsigned levelNew ) const {
    Matrix out( _settings->GetNStates(), _settings->GetNqPEAtRef( levelNew ), 0.0 );
    std::vector<unsigned> qIndex = _settings->GetIndicesQforRef( levelNew );
    unsigned nTotal              = _nTotalForRef[levelOld];

    for( unsigned s = 0; s < _settings->GetNStates(); ++s ) {
        for( unsigned k = 0; k < qIndex.size(); ++k ) {
            for( unsigned i = 0; i < nTotal; ++i ) {
                out( s, k ) += lambda( s, i ) * _phiTildeTrans( i, qIndex[k] );
            }
        }
    }
    return out;
}

void MultiElementSG::EvaluateLambda( Matrix& out, const Matrix& lambda ) const { out = lambda * _phiTildeTrans; }
