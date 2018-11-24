__precompile__
mutable struct Settings
    # grid settings
    # number spatial interfaces
    Nx::Int64;
    # start and end point
    a::Float64;
    b::Float64;
    # grid cell width
    dx::Float64

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
    NXi::Int64;

    # solution values
    uL::Float64;
    uM::Float64;
    uR::Float64;
    x0::Float64;
    x1::Float64;

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

    sigma0::Float64;
    sigma1::Float64;
    lambdaFilter::Float64;
    filterType::String

    # initial condition
    IC;
    solutionExact;
    intervals::Int64;
    closureType

    function Settings()
        Nx = 6000; #2000
        a = 0.0;
        b = 1.0;
        x = linspace(a,b,Nx)
        dx = (b-a)/(Nx-1.0);
        tEnd = 0.01115; # 0.055 0.01115
        N = 2;#5 8
        NXi = 4;
        level = 3; # 3 4
        Nq = 3;#ceil(1.5*N+1);
        uL = 12.0;
        uM = 6.0;
        uR = 1.0;
        x0 = 0.3; # 0.5
        x1 = 0.6;
        cfl = 0.8;
        epsAbs = 0.0;
        sigma0 = 0.2; # 0.2
        sigma1 = 0.1;

        problem = "Burgers"; # Possibilities are Advection and Burgers
        quadratureType = "Gauss"; # Possibilities are Gauss and ClenshawCurtis

        filterType = "noFilter" # Possibilities are noFilter, L1Filter, L2Filter, L1FilterClosure, hSG, hL1Filter, hSGSummed, hLassoSummed
        lambdaFilter = 0.000002; #0.0000001;


        advectionSpeed = 10.0;

        closureType = "BoundedBarrier"  # BoundedBarrier, LogBarrier
        du = 0*1e-4;
        if closureType == "BoundedBarrier"
            uMinus = uR-du;
            uPlus = uL+sigma0+du;
        else
            uMinus = 0.0;
            uPlus = 12.5;
        end

        if problem == "Burgers"
            dt = cfl*dx/uPlus;
        else
            dt = cfl*dx/advectionSpeed;
        end

        epsNewton = 1.0*10.0^(-5);
        NMaxNewton = 100;

        intervals = 10;
        
        # build class
        new(Nx,a,b,dx,tEnd,dt,N,Nq,quadratureType,N,NXi,uL,uM,uR,x0,x1,x,problem,advectionSpeed,uMinus,uPlus,epsNewton,NMaxNewton,epsAbs,sigma0,sigma1,lambdaFilter,filterType,
            (x,xi,eta)->IC3(x,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1),
            (t,x,xi,eta)->IC3Exact(t,x,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1,problem,advectionSpeed),
            intervals,closureType);
    end

end

function IC2(x,xi,sigma,uL,uR,x0,x1)
    y = zeros(size(xi));
    K0 = 12;
    K1 = 1;
    div = 1/(x0^3+3*x0*x1^2-x1^3-3*x1*x0^2);
    a = -2*(K0-K1)*div;
    b = 3*(K0-K1)*(x0+x1)*div;
    c = -6*(K0-K1)*x0*x1*div;
    d = (-K0*x1^3+3*K0*x0*x1^2+K1*x0^3-3*K1*x1*x0^2)*div;
    
    for k = 1:length(xi)
        if x-sigma*xi[k] < x0
            y[k] = K0;
        elseif x-sigma*xi[k] < x1
            y[k] = a*(x-sigma*xi[k])^3+b*(x-sigma*xi[k])^2+c*(x-sigma*xi[k])+d;
        else
            y[k] = K1;
        end
    end
    return y;
end

function IC3Exact(t::Float64,x::Float64,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1,problem,advectionSpeed)
    y = zeros(length(xi));
    for j = 1:length(y);
        tStar = 2.0*(x0-x1)/(uR-uL-sigma0*xi[j]);
        if t < tStar
            v0 = 0.5*(uL+uM+sigma0*xi[j]+sigma1*eta[j]);
            v1 = 0.5*(uR+uM+sigma1*eta[j]);
            x0T = x0+v0*t;
            x1T = x1+v1*t
            if x < x0T
                y[j] = uL+ sigma0*xi[j];
            elseif x < x1T
                y[j] = uM+ sigma1*eta[j];
            else
                y[j] = uR;
            end
        end
        if t >= tStar
            xStar = (x0-x1)*(uM-sigma1*eta+uR)/(uR-uL-sigma0*xi)+x1;
            vStar = 0.5*(uL + sigma0*xi +uR);
            if x < xStar + (t-tStar)*vStar
                y[j] = uL + sigma0*xi[j];
            else
                y[j] = uR;
            end
        end
    end
    return y;
end

function IC1(x,xi,sigma,uL,uR,x0,x1)
    y = zeros(length(xi));
    for j = 1:length(y);
        if x < x0+sigma*xi[j]
            y[j] = uL;
        elseif x < x1+sigma*xi[j]
            y[j] = uL + (uR - uL)*(x-sigma*xi[j]-x0)/(x1-x0);
        else
            y[j] = uR;
        end
    end
    return y;
end

function IC3(x,xi,eta,sigma0,sigma1,uL,uM,uR,x0,x1)
    y = zeros(length(xi));
    for j = 1:length(y);
        if x < x0
            y[j] = uL + sigma0*xi[j];
        elseif x < x1
            y[j] = uM + sigma1*eta[j];
        else
            y[j] = uR;
        end
    end
    return y;
end

function IC1Exact(t::Float64,x::Float64,xi,sigma,uL,uR,x0,x1,problem,advectionSpeed)
    y = zeros(length(xi));
    x0Save = x0;
    x1Save = x1;

    if problem == "Burgers"
        for j = 1:length(y);
            if t >= (x1Save-x0Save)/(uL-uR);
                tS = (x1Save-x0Save)/(uL-uR);
                x0BeforeShock = x0Save+sigma*xi[j] + tS*uL;
                x1BeforeShock = x1Save+sigma*xi[j] + tS*uR;
                x0 = x0BeforeShock + (t-tS)*(uL+uR)*0.5;
                x1 = x0 - 1.0; # ???
            else
                x0 = x0Save+sigma*xi[j] + t*uL;
                x1 = x1Save+sigma*xi[j] + t*uR;
            end
            if x < x0
                y[j] = uL;
            elseif x < x1
                y[j] = uL + (uR - uL)*(x-x0)/(x1-x0);
            else
                y[j] = uR;
            end
        end
    elseif problem == "Advection"
        for j = 1:length(y);
            x0 = x0+sigma*xi[j] + t*advectionSpeed;
            x1 = x1+sigma*xi[j] + t*advectionSpeed;
            if x < x0
                y[j] = uL;
            elseif x < x1
                y[j] = uL + (uR - uL)*(x-x0)/(x1-x0);
            else
                y[j] = uR;
            end
        end
    end
    return y;
end
