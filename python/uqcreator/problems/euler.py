import numpy as np
from problem import Problem

class Euler(Problem):
    def __init__(self, dim, nStates):
        self.dim = dim
        self.nStates = nStates
        self.nQ = 1
        self.gamma = 1.4
    
    def G(self, u, v, nUnit, n, level):        
        y = np.zeros( self.nStates, self.nQ )
        for k in range( self.nQ ):        
            rhoInv = 1.0 / u[0, k]
            uU     = u[1, k] * rhoInv
            vU     = u[2, k] * rhoInv
            p      = ( self.gamma - 1.0 ) * ( u[3, k] - 0.5 * u( 0, k ) * ( uU**2 +  vU**2 ) )
            aU     = np.sqrt( self.gamma * p * rhoInv )

            rhoInv = 1.0 / v( 0, k )
            uV     = v[1, k] * rhoInv
            vV     = v[2, k] * rhoInv
            p      = ( self.gamma - 1.0 ) * ( v( 3, k ) - 0.5 * v( 0, k ) * ( uV**2 + vV**2 ) )
            aV     = np.sqrt( self.gamma * p * rhoInv )

            uUProjected = nUnit[0] * uU + nUnit[1] * vU
            uVProjected = nUnit[0] * uV + nUnit[1] * vV

            lambdaMin = uUProjected - aU
            lambdaMax = uVProjected + aV

            if lambdaMin >= 0:
                y[:,k] = self.F( u[:,k] ) * n
            elif lambdaMax <= 0:
                y[:,k] = self.F( v[:,k] ) * n
            else:
                y[:,k] = ( 1.0 / ( lambdaMax - lambdaMin ) ) * ( ( lambdaMax * self.F( u[:,k] ) - lambdaMin * self.F( v[:,k] ) ) * n + \
                                                                   lambdaMax * lambdaMin * ( v[:,k] - u[:,k] ) * np.norm( n ) )    
        return y

    def F(self, u):
        v1     = u[1] / u[0]
        v2     = u[2] / u[0]
        p      = ( self.gamma - 1.0 ) * ( u[3] - 0.5 * u[0] * ( v1**2 + v2**2 ) )
        
        flux = np.zeros( self.nStates, 2 )
        flux[0,0] = u[1]
        flux[1,0] = u[1] * v1 + p
        flux[2,0] = u[1] * v2
        flux[3,0] = ( u[3] + p ) * v1
        flux[0,1] = u[2]
        flux[1,1] = u[2] * v1
        flux[2,1] = u[2] * v2 + p
        flux[3,1] = ( u[3] + p ) * v2

        return flux

    def computeDt(self, level):
        dtMinTotal = 1e10;
        kEnd = _settings->GetNqPEAtRef( level );

        cfl = _settings->GetCFL();
        for l in range( _settings->GetNMultiElements() ):
            for k in range(kEnd):
                rhoInv = 1.0 / u[0, l, k ]
                uU     = u[ 1, l, k ] * rhoInv
                vU     = u[ 2, l, k ] * rhoInv
                p      = ( self.gamma - 1.0 ) * ( u[ 3, l, k ] - 0.5 * u[ 0, l, k ] * ( uU**2 + vU**2 ) )
                a      = sqrt( self.gamma * p * rhoInv )

                dtMin      = ( cfl / dx ) * min( min( abs( 1.0 / ( vU - a ) ), abs( 1.0 / ( vU + a ) ) ), \
                                                 min( abs( 1.0 / ( uU + a ) ), abs( 1.0 / ( uU - a ) ) ) )
                dtMinTotal = min( dtMin, dtMinTotal )
            }
        }

        return dtMinTotal;

    def IC(self, x, xi):
        gamma       = 1.4;
        R           = 287.87;
        T           = 273.15;
        Ma          = 0.8;
        AoA         = 1.25;
        AoAScaling  = 1.0;
        p           = 101325.0;
        rhoFarfield = p / ( R * T );        

        if( len(xi) > 1 ) {
            Ma = Ma - sigma[1];
            Ma = Ma + xi[1] * sigma[1];
        }

        a = sqrt( gamma * R * T );
        uMax  = Ma * a;
        angle = ( AoA + AoAScaling * _sigma[0] * xi[0] ) * ( 2.0 * M_PI ) / 360.0;
        uF    = uMax * cos( angle );
        vF    = uMax * sin( angle );

        y = np.zeros(self.nStates)
        y[0]                  = rhoFarfield;
        y[1]                  = rhoFarfield * uF;
        y[2]                  = rhoFarfield * vF;
        kineticEnergyL = 0.5 * rhoFarfield * ( pow( uF, 2 ) + pow( vF, 2 ) );
        innerEnergyL   = ( p / ( rhoFarfield * ( gamma - 1 ) ) ) * rhoFarfield;
        y[3]                  = kineticEnergyL + innerEnergyL;
        return y;       
        