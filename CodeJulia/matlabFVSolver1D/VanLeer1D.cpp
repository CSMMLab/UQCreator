#include "VanLeer1D.h"

VanLeer1D::VanLeer1D( PhysicalProblem* physicalProblem ) : Limiter( physicalProblem ) {}

VanLeer1D::~VanLeer1D() {}

void VanLeer1D::CalculateLimitedState1D( int index, Mesh* mesh, Vector* limitLeft, Vector* limitRight ) {
    double theta, phi, sigma;
    double dxInv = 1 / mesh->GetDx();
    double eps   = 1e-8;
    Vector *left, *current, *right;
    current = mesh->GetCell( index )->GetState();
    left    = mesh->GetLeftCell( index )->GetState();
    right   = mesh->GetRightCell( index )->GetState();
    for( unsigned int i = 0; i < _physicalProblem->GetStateDim(); i++ ) {
        theta              = ( ( *current )[i] - ( *left )[i] ) / ( ( *right )[i] - ( *current )[i] + eps );
        phi                = ( theta + std::fabs( theta ) ) / ( 1 + std::fabs( theta ) );
        sigma              = ( ( *right )[i] - ( *current )[i] ) * dxInv * phi;
        ( *limitLeft )[i]  = ( *current )[i] - 0.5 * mesh->GetDx() * sigma;
        ( *limitRight )[i] = ( *current )[i] + 0.5 * mesh->GetDx() * sigma;
    }
}

void VanLeer1D::t1_CalculateLimitedState1D(
    int index, Mesh* mesh, Mesh* t1_mesh, Vector* limitLeft, Vector* t1_limitLeft, Vector* limitRight, Vector* t1_limitRight ) {
    double t1_theta, t1_phi, t1_sigma;
    double theta, phi, sigma;
    double dxInv = 1 / mesh->GetDx();
    double eps   = 1e-8;
    Vector *t1_left, *t1_current, *t1_right;
    Vector *left, *current, *right;

    t1_current = t1_mesh->GetCell( index )->GetState();
    current    = mesh->GetCell( index )->GetState();

    t1_left = t1_mesh->GetLeftCell( index )->GetState();
    left    = mesh->GetLeftCell( index )->GetState();

    t1_right = t1_mesh->GetRightCell( index )->GetState();
    right    = mesh->GetRightCell( index )->GetState();

    for( unsigned int i = 0; i < _physicalProblem->GetStateDim(); i++ ) {

        t1_theta = ( ( *right )[i] - ( *left )[i] + eps ) / pow( ( ( *right )[i] - ( *current )[i] + eps ), 2 ) * ( *t1_current )[i] -
                   ( ( *current )[i] - ( *left )[i] + eps ) / pow( ( ( *right )[i] - ( *current )[i] + eps ), 2 ) * ( *t1_right )[i] -
                   ( 1 / ( ( *right )[i] - ( *current )[i] + eps ) ) * ( *t1_left )[i];
        theta = ( ( *current )[i] - ( *left )[i] ) / ( ( *right )[i] - ( *current )[i] + eps );

        if( theta < 0 ) {
            t1_phi = 0;
        }
        else {
            t1_phi = 2 / pow( ( theta + 1 ), 2 ) * t1_theta;
        }
        phi = ( theta + std::fabs( theta ) ) / ( 1 + std::fabs( theta ) );

        t1_sigma = phi * dxInv * ( *t1_right )[i] - phi * dxInv * ( *t1_current )[i] + ( ( *right )[i] - ( *current )[i] ) * dxInv * t1_phi;
        sigma    = ( ( *right )[i] - ( *current )[i] ) * dxInv * phi;

        ( *t1_limitLeft )[i] = ( *t1_current )[i] - 0.5 * mesh->GetDx() * t1_sigma;
        ( *limitLeft )[i]    = ( *current )[i] - 0.5 * mesh->GetDx() * sigma;

        ( *t1_limitRight )[i] = ( *t1_current )[i] + 0.5 * mesh->GetDx() * t1_sigma;
        ( *limitRight )[i]    = ( *current )[i] + 0.5 * mesh->GetDx() * sigma;
    }
}
