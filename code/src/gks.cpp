void get_conserved(double *w, double prim[3], double gam)
{
    w[0] = prim[0];
    w[1] = prim[0] * prim[1];
    w[2] = 0.5 * prim[0] / prim[2] / (gam - 1.0) + 0.5 * prim[0] * prim[1] * prim[1];
}

void get_primitive(double *prim, double w[3], double gam)
{
    prim[0] = w[0];
    prim[1] = w[1] / w[0];
    prim[2] = 0.5 * w[0] / (gam - 1.0) / (w[2] - 0.5 * w[1] * w[1] / w[0]);
}

double get_tau(double density_left, double density_right, double density0, 
               double lambda_left, double lambda_right, double lambda0, 
               double mu, double dt)
{
    double tau_num = 1.0 * fabs(density_left / lambda_left - density_right / lambda_right) 
                     / fabs(density_left / lambda_left + density_right / lambda_right) * dt;
    
    double tau_ns = 2.0 * mu * lambda0 / density0;
    
    return tau_num + tau_ns;
    
}

void calc_moment(double *Mu, double *Mu_L, double *Mu_R, double *Mxi, double prim[3], double inK)
{
    
    const double PI = 3.1415926535;

    // moments of normal velocity
    Mu_L[0] = 0.5 * erfc(-sqrt(prim[2]) * prim[1]);
    Mu_L[1] = prim[1] * Mu_L[0] + 0.5 * exp(-prim[2] * pow(prim[1], 2)) / sqrt(PI * prim[2]);
    Mu_R[0] = 0.5 * erfc(sqrt(prim[2]) * prim[1]);
    Mu_R[1] = prim[1] * Mu_R[0] - 0.5 * exp(-prim[2] * pow(prim[1], 2)) / sqrt(PI * prim[2]);

    for (int i=2;i<=6;i++)
    {
        Mu_L[i] = prim[1] * Mu_L[i-1] + 0.5 * (i - 1) * Mu_L[i-2] / prim[2];
        Mu_R[i] = prim[1] * Mu_R[i-1] + 0.5 * (i - 1) * Mu_R[i-2] / prim[2];
    }

    for (int i=0;i<=6;i++)
    {
        Mu[i] = Mu_L[i] + Mu_R[i];
    }

    // moments of inner degrees of freedom
    Mxi[0] = 1.0;
    Mxi[1] = 0.5 * inK / prim[2];
    Mxi[2] = (pow(inK, 2) + 2.0 * inK) / (4.0 * pow(prim[2], 2));
    
}

void moment_uv(double *moment_uv, double Mu[7], double Mxi[3], int alpha, int delta)
{
    moment_uv[0] = Mu[alpha] * Mxi[delta/2];
    moment_uv[1] = Mu[alpha+1] * Mxi[delta/2];
    moment_uv[2] = 0.5 * (Mu[alpha+2] * Mxi[delta/2] + Mu[alpha] * Mxi[(delta+2)/2]);
}

void moment_au(double *moment_au, double a[3], double Mu[7], double Mxi[3], int alpha)
{
    double t0[3], t1[3], t2[3], t3[3];
    
    moment_uv(t0, Mu, Mxi, alpha+0, 0);
    moment_uv(t1, Mu, Mxi, alpha+1, 0);
    moment_uv(t2, Mu, Mxi, alpha+2, 0);
    moment_uv(t3, Mu, Mxi, alpha+0, 2);
    
    for (int i=0;i<=2;i++)
    {
        moment_au[i] = a[0] * t0[i] + a[1] * t1[i] + 0.5 * a[2] * t2[i] + 0.5 * a[2] * t3[i];
    }
}

void GKS( double *flux, double u[3], double v[3], double gam, double mu, double dt ) 
{    
    // gas property
    double inK = (3.0 - gam) / (gam - 1.0);
    
    // left interface
    double wL[3], primL[3];
    
    wL[0] = u[0];
    wL[1] = u[1];
    wL[2] = u[2];
    
    get_primitive(primL, wL, gam);
    
    double pL = ( gam - 1.0 ) * ( wL[2] - 0.5 * wL[0] * pow( wL[1], 2 ) );
    double ssL = sqrt( gam * pL / primL[0] );
    
    // right interface
    double wR[3], primR[3];
    
    wR[0] = v[0];
    wR[1] = v[1];
    wR[2] = v[2];
    
    get_primitive(primR, wR, gam);
    
    double pR = ( gam - 1.0 ) * ( wR[2] - 0.5 * wR[0] * pow( wR[1], 2 ) );
    double ssR = sqrt( gam * pR / primR[0] );
    
    // central interface
    double Mu[7], MuL[7], MuR[7], Mxi[3];
    double Mau[3], MauL[3], MauR[3];
    
    calc_moment(Mu, MuL, MuR, Mxi, primL, inK);
    moment_uv(MauL, Mu, Mxi, 0, 0);
    
    calc_moment(Mu, MuL, MuR, Mxi, primR, inK);
    moment_uv(MauR, Mu, Mxi, 0, 0);
    
    double w[3], prim[3], tau;
    for (int i=0;i<=2;i++)
    {
        w[i] = primL[0] * MauL[i] + primR[1] * MauR[i];
    }
    
    get_primitive(prim, w, gam);
    
    tau = get_tau(primL[0], primR[0], prim[0], primL[2], primR[2], prim[2], mu, dt);
    
    // time integral terms
    double Mt[5];
    
    Mt[3] = tau * (1.0 - exp(-dt / tau));
    Mt[4] = -tau * dt * exp(-dt / tau) + tau * Mt[3];
    Mt[0] = dt - Mt[3];
    Mt[1] = -tau * Mt[0] + Mt[4];
    Mt[2] = dt * dt / 2.0 - tau * Mt[0];
    
    // calculate the flux of conservative variables related to g0
    calc_moment(Mu, MuL, MuR, Mxi, prim, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    
    for (int i=0;i<=2;i++)
    {
        flux[i] = Mt[0] * primL[0] * Mau[i];
    }

    // calculate the flux of conservative variables related to f0
    calc_moment(Mu, MuL, MuR, Mxi, primL, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    for (int i=0;i<=2;i++)
    {
        flux[i] += Mt[3] * primL[0] * Mau[i];
    }
    
    calc_moment(Mu, MuL, MuR, Mxi, primR, inK);
    moment_uv(Mau, Mu, Mxi, 1, 0);
    for (int i=0;i<=2;i++)
    {
        flux[i] += Mt[3] * primR[0] * Mau[i];
    }
    
    // final flux
    for (int i=0;i<=2;i++)
    {
        flux[i] = flux[i] / dt;
    }

}
