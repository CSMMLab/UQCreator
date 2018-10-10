__precompile__
include("analytic_sod.jl")
mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    Ny::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    c::Float64;
    d::Float64;
    # grid cell width
    dx::Float64
    dy::Float64

    # time settings
    # end time
    tEnd::Float64;
    # time increment
    dt::Float64;
    # Runge Kutta stages of time integrator
    rkStages::Int64;

    # number of quadrature points
    Nq::Int64;
    # definition of quadrature type
    quadratureType::String;

    # maximal polynomial degree
    N::Int64;

    # solution values
    rhoL::Float64;
    rhoR::Float64;
    pL::Float64;
    pR::Float64;
    uL::Float64;
    uR::Float64;

    # initial shock position at xi = 0
    x0::Float64;

    # grid
    x

    # problem definitions
    problem::String;

    # advection speed if advection is defined as problem
    advectionSpeed::Float64;

    # imposed solution bounds
    uMinus::Float64;
    uPlus::Float64;

    # settings for newton solver
    epsNewton::Float64;
    NMaxNewton::Int64;

    # settings for QBC abs approximation
    epsAbs::Float64;

    sigma::Float64;
    lambdaFilter::Float64;
    filterType::String

    # initial condition
    IC;
    solutionExact;

    # constants for Euler Problem
    gamma::Float64;
    SpecificR::Float64;
    states::Int;
    dimensions::Int;

    function Settings()
        Nx = 50; #400
        Ny = 50; #400
        a = 0.0;
        b = 1.0;
        c = 0.0;
        d = 1.0;
        x = linspace(a,b,Nx+1)
        y = linspace(c,d,Ny+1)
        dx = x[2]-x[1];
        dy = y[2]-y[1];
        N = 5;#9;
        Nq = ceil(1.5*N+0.5);

        x0 = 0.5; # 0.5
        cfl = 0.5;
        epsAbs = 0.0;

        problem = "Burgers"; # Possibilities are Advection and Burgers
        quadratureType = "Gauss"; # Possibilities are Gauss and ClenshawCurtis

        filterType = "L2Filter" # Possibilities are noFilter, L1Filter, L2Filter, L1FilterClosure
        lambdaFilter = 0.000003;


        advectionSpeed = 10.0;

        dt = 0.2/max(Nx,Ny);
        tEnd = 0.8; # 0.35

        du = 0.001;

        epsNewton = 5.0*10.0^(-5);
        NMaxNewton = 1000;

        sigma = 0.1; # 0.05 0.2
        factor = 3;
        gamma = 1.4;
        SpecificR = 462.0; # only for temperature calculation
        rhoL = 1.0; rhoR = 0.8; #0.3;
        pL = 1.0; pR = 0.3;
        uL = 0.0; uR = 0.0;
        vL = 0.0; vR = 0.0;

        uMinus = uR;
        uPlus = uL;

        states = 4; # number of conserved states i.e. density, momentum, energy
        dimensions = 2; # dimension of physical space
        
        # build class
        new(Nx,Ny,a,b,c,d,dx,dy,tEnd,dt,N,Nq,quadratureType,N,rhoL,rhoR,pL,pR,uL,uR,x0,x,problem,advectionSpeed,uMinus,uPlus,epsNewton,NMaxNewton,epsAbs,sigma,lambdaFilter,filterType,
            (s,x,xi)->ICEuler(s,x,xi,sigma,rhoL,rhoR,pL,pR,uL,uR,vL,vR,x0,gamma),
            (t,xV,xi)->analytic_sod(t,x0+sigma*xi,rhoL,pL,uL,rhoR,pR,uR,gamma,xV),
            gamma,SpecificR,states,dimensions);
    end

end

function ICEuler(s,x::Array{Float64,1},xi,sigma,rhoL,rhoR,pL,pR,uL,uR,vL,vR,x0,gamma)
    y = zeros(size(xi));
    indexX = 1;
    if s == 1        
        for j = 1:length(y);
            if x[indexX] < x0+sigma*xi[j];
                y[j] = rhoL;
            else
                y[j] = rhoR;
            end
        end
    elseif s == 2
        rhouL=rhoL*uL;
        rhouR=rhoR*uR;
        for j = 1:length(y);
            if x[indexX] < x0+sigma*xi[j];
                y[j] = rhouL;
            else
                y[j] = rhouR;
            end
        end
    elseif s == 3
        rhovL=rhoL*vL;
        rhovR=rhoR*vR;
        for j = 1:length(y);
            if x[indexX] < x0+sigma*xi[j];
                y[j] = rhovL;
            else
                y[j] = rhovR;
            end
        end
    elseif s == 4
        kineticEnergyL=0.5*rhoL*(uL^2+vR^2);
        innerEnergyL=(pL/(rhoL*(gamma-1)))*rhoL;
        rhoeL=kineticEnergyL+innerEnergyL;
        kineticEnergyR=0.5*rhoR*(uR^2+vR^2);
        innerEnergyR=(pR/(rhoR*(gamma-1)))*rhoR;
        rhoeR=kineticEnergyR+innerEnergyR;
        for j = 1:length(y);
            if x[indexX] < x0+sigma*xi[j];
                y[j] = rhoeL;
            else
                y[j] = rhoeR;
            end
        end
    end

    return y;
end

function ICShock(s,x::Array{Float64,1},xi,sigma,rhoL,rhoR,pL,pR,uL,uR,vL,vR,x0,gamma)
    y = zeros(size(xi));
    width = 0.05;
    xM = [x0;x0];
    if s == 1        
        for j = 1:length(y);
            if norm(x-xM) < sigma*xi[j]+width
            #if x[1] > x0+sigma*xi[j]-width && x[1] < x0+sigma*xi[j]+width && x[2] > x0+sigma*xi[j]-width && x[2] < x0+sigma*xi[j]+width;
                    y[j] = rhoL;
            else
                y[j] = rhoR;
            end
        end
    elseif s == 2
        rhouL=rhoL*uL;
        rhouR=rhoR*uR;
        for j = 1:length(y);
            if norm(x-xM) < sigma*xi[j]+width
                y[j] = rhouL;
            else
                y[j] = rhouR;
            end
        end
    elseif s == 3
        rhovL=rhoL*vL;
        rhovR=rhoR*vR;
        for j = 1:length(y);
            if norm(x-xM) < sigma*xi[j]+width
                y[j] = rhovL;
            else
                y[j] = rhovR;
            end
        end
    elseif s == 4
        kineticEnergyL=0.5*rhoL*(uL^2+vR^2);
        innerEnergyL=(pL/(rhoL*(gamma-1)))*rhoL;
        rhoeL=kineticEnergyL+innerEnergyL;
        kineticEnergyR=0.5*rhoR*(uR^2+vR^2);
        innerEnergyR=(pR/(rhoR*(gamma-1)))*rhoR;
        rhoeR=kineticEnergyR+innerEnergyR;
        for j = 1:length(y);
            if norm(x-xM) < sigma*xi[j]+width
                y[j] = rhoeL;
            else
                y[j] = rhoeR;
            end
        end
    end

    return y;
end
