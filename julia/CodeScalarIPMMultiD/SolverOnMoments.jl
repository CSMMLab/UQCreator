__precompile__
include("quadrature.jl")
include("Basis.jl")

using ProgressMeter
using TensorOperations

struct SolverOnMoments
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # quadrature
    q::Quadrature;

    # SolverOnMoments settings
    settings::Settings;

    # spatial basis functions
    basis::Basis;

    # filter strength
    lambdaFilter::Float64;

    # Moment flux matrix
    B::Array{Float64,3}

    # constructor
    function SolverOnMoments(settings)
        settings.Nq = 100;
        x = linspace(settings.a,settings.b,settings.Nx)
        q = Quadrature(settings.Nq,"Gauss");
        basis = Basis(q,settings);

        N = settings.N;
        B = zeros(N,N,N);
        for i = 1:N
            for n = 1:N
                for m = 1:N
                    B[i,n,m] = IntegralVec(q,basis.PhiQuad[:,i].*basis.PhiQuad[:,n].*basis.PhiQuad[:,m]*0.5,-1.0,1.0);
                end
            end
        end

        new(x,q,settings,basis,s.lambdaFilter,B);
    end
end

function numFlux(obj::SolverOnMoments,u::Float64,v::Float64)
    return 0.5*(fBurgers(u)+fBurgers(v)-(obj.settings.dx/obj.settings.dt).*(v-u));
end

# Local Lax-Friedrichs (LLF) flux for Burgers
function numFluxMoments(obj::SolverOnMoments,u::Array{Float64,1},v::Array{Float64,1})
    y = zeros(obj.settings.N);
    B = obj.B;
    #for i = 1:obj.settings.N
    #    for n = 1:obj.settings.N
    #        for m = 1:obj.settings.N
    #            y[i] = y[i] + v[m]*obj.B[m,n,i]*v[n] + u[m]*obj.B[m,n,i]*u[n];
    #        end
    #    end
    #end

    @tensoropt y[i] = v[m]*B[m,n,i]*v[n]+  u[m]*B[m,n,i]*u[n];

    #for i = 1:obj.settings.N
        #y[i] = y[i] + sum(v.*(obj.B[:,:,i]*v)) + sum(u.*(obj.B[:,:,i]*u));
        #y[i] = y[i] + dot(v,obj.B[:,:,i]*v) + dot(u,obj.B[:,:,i]*u)
        #y[i] =  dot(v,obj.B[:,:,i]*v) + dot(u,obj.B[:,:,i]*u)
        #sum(x.*y)
    #end

    #y =  [v'*obj.B[:,:,i]*v+u'*obj.B[:,:,i]*u for i = 1:obj.settings.N]

    return 0.25*(y-(obj.settings.dx/obj.settings.dt)*(v-u));
end

function SetupIC(obj::SolverOnMoments)
    u = zeros(obj.settings.N,obj.settings.Nx);
    for j = 1:obj.settings.Nx
        for i = 1:obj.settings.N
            uVals = obj.settings.IC(obj.x[j],obj.q.xi);
            u[i,j] = IntegralVec(obj.q, uVals.*obj.basis.PhiQuad[:,i]*0.5,-1.0,1.0);
        end
    end
    return u;
end

function Solve(obj::SolverOnMoments)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    N = obj.settings.N; # number of moments

    # Set up initial condition
    u = SetupIC(obj);
    for j = 1:Nx
        u[:,j] = Filter(obj,u[:,j]);
    end
    uNew = zeros(size(u));
    
    Nt = round(tEnd/dt);
    
    # time loop
    @showprogress 0.1 "Progress " for n = 1:Nt
        #println("t = ",t);

        # Update time by dt
        for j = 2:(Nx-1)
            uNew[:,j] = u[:,j]-dt/dx*(numFluxMoments(obj,u[:,j],u[:,j+1])-numFluxMoments(obj,u[:,j-1],u[:,j]));
        end

        for j = 2:(Nx-1)
            u[:,j] = Filter(obj,uNew[:,j]);
        end
        t = t+dt;
    end

    # return end time and solution
    return t, u;

end

function Filter(obj::SolverOnMoments,u::Array{Float64,1})
    if obj.settings.filterType == "L2Filter"
        for i = 1:obj.settings.N
            u[i] = u[i]/(1.0+obj.lambdaFilter*i^2*(i-1)^2);
        end
    elseif obj.settings.filterType == "L1Filter"
        for i = 1:obj.settings.N
            scL1 = 1.0-4.0*obj.lambdaFilter*i*(i-1)*obj.basis.PhiL1[i]/abs(u[i]); # 4.0 term to scale filter strength, can be left out
            if scL1 < 0
                scL1 = 0;
            end
            u[i] = scL1*u[i];
        end
    elseif obj.settings.filterType == "L1FilterClosure"
        N = obj.settings.N;
        lambdaFilter = abs(u[N])/(N*(N-1)*obj.basis.PhiL1[N]);
        for i = 2:obj.settings.N
            scL1 = 1.0-lambdaFilter*i*(i-1)*obj.basis.PhiL1[i]/abs(u[i]);
            if scL1 < 0 || abs(u[i]) < 1.0e-7
                scL1 = 0;
            end
            u[i] = scL1*u[i];
        end
    end
    return u;
end

function f(obj::SolverOnMoments,u)
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