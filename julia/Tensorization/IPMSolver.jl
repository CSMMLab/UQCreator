__precompile__
include("quadrature.jl")
include("Basis.jl")
include("closure.jl")

using ProgressMeter

struct IPMSolver
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # quadrature
    q::Quadrature;

    # IPMSolver settings
    settings::Settings;

    # spatial basis functions
    basis::Basis;

    # filter strength
    lambdaFilter::Float64;

    # IPM Closure
    closure::Closure

    # constructor
    function IPMSolver(settings)
        x = linspace(settings.a,settings.b,settings.Nx)
        q = Quadrature(settings.Nq,"Gauss");
        basis = Basis(q,settings);
        closure = Closure(settings,basis,q)

        new(x,q,settings,basis,s.lambdaFilter,closure);
    end
end

# Local Lax-Friedrichs (LLF) flux for Burgers
function numFlux(obj::IPMSolver,u::Float64,v::Float64)
    beta = max(abs(u),abs(v)); # max of abs(f'(u)) 
    return 0.5*(fBurgers(u)+fBurgers(v)-beta.*(v-u));
end

function SetupIC(obj::IPMSolver)
    u = zeros(obj.settings.N,obj.settings.Nx);
    for j = 1:obj.settings.Nx
        for i = 1:obj.settings.N
            uVals = obj.settings.IC(obj.x[j],obj.q.xi);
            u[i,j] = IntegralVec(obj.q, uVals.*obj.basis.PhiQuad[:,i]*0.5,-1.0,1.0);
        end
    end
    return u;
end

function Solve(obj::IPMSolver)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    N = obj.settings.N; # number of moments

    # Set up initial condition
    u = SetupIC(obj);
    v = zeros(size(u));
    vNew = zeros(size(u));
    for j = 1:Nx
        vNew[:,j] = Solve(obj.closure,u[:,j],v[:,j]);
        v[:,j] = vNew[:,j];
        u[:,j] = ComputeMoments(obj.closure,vNew[:,j]);
    end
    uNew = zeros(size(u));
    
    Nt = round(tEnd/dt);
    
    # time loop
    @showprogress 0.1 "Progress " for n = 1:Nt
        #println("t = ",t);

        # Update time by dt
        for j = 2:(Nx-1)
            numFRight = numFlux.(obj, Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j])),Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j+1])));
            numFLeft = numFlux.(obj, Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j-1])),Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j])));
            for i = 1:N
                uNew[i,j] = u[i,j]-dt/dx*IntegralVec(obj.q,(numFRight-numFLeft).*obj.basis.PhiQuad[:,i]*0.5,-1.0,1.0);
            end
        end

        for j = 2:(Nx-1)
            vNew[:,j] = Solve(obj.closure,uNew[:,j],v[:,j]);
            v[:,j] = vNew[:,j];
            u[:,j] = ComputeMoments(obj.closure,vNew[:,j]);
        end
        t = t+dt;
    end

    # return end time and solution
    return t, u, v;

end

function f(obj::IPMSolver,u)
    if obj.settings.problem == "Advection"
        return obj.settings.advectionSpeed*u;
    else
        return 0.5*u.^2.0;
    end
end

# physical flux Scalar valued Burgers
function fBurgers(u::Float64)
    return 0.5*u^2.0;
end

# physical flux Vector valued Burgers
function fBurgers(u::Array{Float64,1})
    return 0.5*u.^2.0;
end

# physical flux Scalar valued Advection
function fAdvection(advectionSpeed::Float64, u::Float64)
    return advectionSpeed*u;
end

# physical flux Vector valued Advection
function fAdvection(advectionSpeed::Float64, u::Array{Float64,1})
    return advectionSpeed*u;
end