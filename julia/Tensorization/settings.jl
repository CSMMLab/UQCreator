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

    # solution values
    uL::Float64;
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

    sigma::Float64;
    lambdaFilter::Float64;
    filterType::String
    barrierRefine::Float64;

    # initial condition
    IC;
    solutionExact;

    function Settings()
        Nx = 1000; #200
        a = 0.0;
        b = 3.0;
        x = linspace(a,b,Nx)
        dx = (b-a)/(Nx-1.0);
        tEnd = 0.14; # 0.11
        N = 2;#15;
        Nq = 3;
        uL = 12.0;
        uR = 3.0;
        x0 = 0.5; # 0.5
        x1 = 1.5;
        cfl = 0.5;
        epsAbs = 0.0;

        problem = "Burgers"; # Possibilities are Advection and Burgers
        quadratureType = "Gauss"; # Possibilities are Gauss and ClenshawCurtis

        filterType = "noFilter" # Possibilities are noFilter, L1Filter, L2Filter, L1FilterClosure
        lambdaFilter = 0.000001;

        advectionSpeed = 10.0;

        du = 0.001;
        uMinus = uR-du;
        uPlus = uL+du;

        if problem == "Burgers"
            dt = cfl*dx/uL;
        else
            dt = cfl*dx/advectionSpeed;
        end


        epsNewton =5e-5; #1.0*dx^2;#
        NMaxNewton = 1000;

        sigma = 0.2; # 0.2

        barrierRefine = 1e-5;
        
        
        # build class
        new(Nx,a,b,dx,tEnd,dt,N,Nq,quadratureType,N,uL,uR,x0,x1,x,problem,advectionSpeed,
            uMinus,uPlus,epsNewton,NMaxNewton,epsAbs,sigma,lambdaFilter,filterType,barrierRefine,
            (x,xi)->IC1(x,xi,sigma,uL,uR,x0,x1),
            (t,x,xi)->IC1Exact(t,x,xi,sigma,uL,uR,x0,x1,problem,advectionSpeed));
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

function IC2Exact(t::Float64,x::Float64,xi,sigma,uL,uR,x0,x1,problem,advectionSpeed)
    y = zeros(size(xi));
    K0 = 12;
    K1 = 1;
    div = 1/(x0^3+3*x0*x1^2-x1^3-3*x1*x0^2);
    a = -2*(K0-K1)*div;
    b = 3*(K0-K1)*(x0+x1)*div;
    c = -6*(K0-K1)*x0*x1*div;
    d = (-K0*x1^3+3*K0*x0*x1^2+K1*x0^3-3*K1*x1*x0^2)*div;
    
    for k = 1:length(xi)
        if x-sigma*xi[k] < x0+t*advectionSpeed
            y[k] = K0;
        elseif x-sigma*xi[k] < x1+t*advectionSpeed
            y[k] = a*(x-sigma*xi[k]-t*advectionSpeed)^3+b*(x-sigma*xi[k]-t*advectionSpeed)^2+c*(x-sigma*xi[k]-t*advectionSpeed)+d;
        else
            y[k] = K1;
        end
    end
    return y;
end

function IC1(x,xi,sigma,uL,uR,x0,x1)
    y = zeros(size(xi));
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
