#include "eulerclosure2d.h"

EulerClosure2D::EulerClosure2D( Problem* problem ) : Closure( problem ), _gamma( problem->GetGamma() ) {}

EulerClosure2D::~EulerClosure2D() {}

void EulerClosure2D::U( Vector& out, const Vector& Lambda ) {
    double v1 = Lambda[0]; double v2 = Lambda[1]; double v3 = Lambda[2]; double v4 = Lambda[3];
    double expTerm = (-exp(((pow(v2,2) + pow(v3,2) - 2.0*v1*v4 - 2.0*v4*_gamma)/(2.0*v4)))*v4)^(1.0/(1.0 + _gamma));
    out[0] = expTerm;
    out[1] = -((v2*expTerm)/v4);
    out[2] = -((v3*expTerm)/v4);
    out[3] = -((expTerm*(-pow(v2,2) - pow(v3,2) + 2.0*v4))/(2.0*pow(v4,2)));
}

Matrix EulerClosure2D::U( const Matrix& Lambda ) {
    double expTerm, v1, v2, v3, v4;
    Matrix y( _nStates, Lambda.columns(), 0.0 );
    for( unsigned k = 0; k < Lambda.columns(); ++k ) {
        v1 = Lambda(0,k); v2 = Lambda(1,k); v3 = Lambda(2,k); v4 = Lambda(3,k);
        expTerm = (-exp(((pow(v2,2) + pow(v3,2) - 2.0*v1*v4 - 2.0*v4*_gamma)/(2.0*v4)))*v4)^(1.0/(1.0 + _gamma));
        y(0,k) = expTerm;
        y(1,k) = -((v2*expTerm)/v4);
        y(2,k) = -((v3*expTerm)/v4);
        y(3,k) = -((expTerm*(-pow(v2,2) - pow(v3,2) + 2.0*v4))/(2.0*pow(v4,2)));
    }

    return y;
}

void EulerClosure2D::DU( Matrix& y, const Vector& Lambda ) {
    double v1 = Lambda[0]; double v2 = Lambda[1]; double v3 = Lambda[2]; double v4 = Lambda[3];
    double expTerm = pow(-exp(((pow(v2,2) + pow(v3,2) - 2.0*v4*(v1 + _gamma))/(2.0*v4)))*v4,(1.0/(1.0 + _gamma)))*1.0/(1.0 + _gamma);
    double vTerm = pow(v2,2) + pow(v3,2) + 2.0*v4*_gamma;
    y(0,0) = -(expTerm);
    y(0,1) = (v2*expTerm)/(v4);
    y(0,2) = (v3*expTerm)/(v4);
    y(0,3) = -(((pow(v2,2) + pow(v3,2) - 2.0*v4)*expTerm)/(2.0*pow(v4,2)));
    y(1,0) = (v2*expTerm)/(v4);
    y(1,1) =  -((expTerm*(pow(v2,2) + v4 + v4*_gamma))/(pow(v4,2)));
    y(1,2) =  -((v2*v3*expTerm)/(pow(v4,2)));
    y(1,3) =  (v2*expTerm*vTerm)/(2*pow(v4,3));
    y(2,0) =  (v3*expTerm)/(v4);
    y(2,1) =  -((v2*v3*expTerm)/(pow(v4,2)));
    y(2,2) =  -((expTerm*(pow(v3,2) + v4 + v4*_gamma))/(pow(v4,2)));
    y(2,3) =  (v3*expTerm *vTerm)/(2.0*pow(v4,3) );
    y(3,0) =  -(((pow(v2,2) + pow(v3,2) - 2.0*v4)*expTerm)/(2.0*pow(v4,2)));
    y(3,1) =  (v2*expTerm *vTerm)/(2.0*pow(v4,3) );
    y(3,2) =  (v3*expTerm *vTerm)/(2.0*pow(v4,3));
    y(3,3) =  -(expTerm*(pow(v2,4) + pow(v3,4) + 4.0*pow(v3,2)*v4*_gamma - 4.0*pow(v4,2)*_gamma + 2.0*pow(v2,2)*(pow(v3,2) + 2.0*v4*_gamma)))/(4.0*pow(v4,4));
}
