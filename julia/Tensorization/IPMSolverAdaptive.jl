__precompile__
include("quadrature.jl")
include("Basis.jl")
include("closure.jl")

using ProgressMeter

struct IPMSolverAdaptive
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # quadrature
    q::Quadrature;

    # quadrature fine
    qF::Quadrature;

    # IPMSolver settings
    settings::Settings;

    # spatial basis functions
    basis::Basis;
 
    # spatial basis functions fine
    basisF::Basis;

    # IPM Closure
    closure::Closure

    # IPM Closure fine
    closureF::Closure

    # refine Array
    refine::Array{Int,2}

    # Lagrange interpolation matrix
    L::Array{Float64,2}


    # constructor
    function IPMSolverAdaptive(settings)
        x = linspace(settings.a,settings.b,settings.Nx)
        q = Quadrature(settings.N,"Gauss");
        basis = Basis(q,settings);
        closure = Closure(settings,basis,q);
        qF = Quadrature(settings.Nq,"Gauss");
        basisF = Basis(qF,settings);
        closureF = Closure(settings,basisF,qF);

        Nt = round(settings.tEnd/settings.dt);
        refine = zeros(Nt+1,settings.Nx);

        L = ones(settings.Nq,settings.N);
        for k = 1:settings.Nq
            for i = 1:settings.N
                for l = 1:settings.N
                    if l != i
                        L[k,i] = L[k,i]*(q.xi[l]-qF.xi[k])/(q.xi[l]-q.xi[i]);
                    end
                end
            end
        end

        new(x,q,qF,settings,basis,basisF,closure,closureF,refine,L);
    end
end

# Local Lax-Friedrichs (LLF) flux for Burgers
function numFlux(obj::IPMSolverAdaptive,u::Float64,v::Float64)
    beta = max(abs(u),abs(v)); # max of abs(f'(u)) 
    return 0.5*(fBurgers(u)+fBurgers(v)-beta.*(v-u));
end

function SetupIC(obj::IPMSolverAdaptive)
    u = zeros(obj.settings.N,obj.settings.Nx);
    uQ = zeros(obj.settings.N,obj.settings.Nx);
    for j = 1:obj.settings.Nx
        uQ[:,j] = obj.settings.IC(obj.x[j],obj.q.xi);
        for i = 1:obj.settings.N
            uVals = obj.settings.IC(obj.x[j],obj.qF.xi);
            u[i,j] = IntegralVec(obj.qF, uVals.*obj.basisF.PhiQuad[:,i]*0.5,-1.0,1.0);
        end
    end
    return u, uQ;
end

function Solve(obj::IPMSolverAdaptive)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    N = obj.settings.N; # number of moments

    # Set up initial condition
    u, uQ = SetupIC(obj);
    v = zeros(size(u));
    vNew = zeros(size(u));
    for j = 1:Nx
        vNew[:,j] = Solve(obj.closureF,u[:,j],v[:,j]);
        v[:,j] = vNew[:,j];
        u[:,j] = ComputeMoments(obj.closureF,vNew[:,j]);
    end
    uNew = deepcopy(u);
    uQNew = deepcopy(uQ);

    UpdateRefinement(obj,u,1);
    
    Nt = round(tEnd/dt);

    # time loop
    @showprogress 0.1 "Progress " for n::Int = 1:Nt

        # Update time by dt
        for j = 2:(Nx-1)
            if obj.refine[n,j] == 0 || obj.refine[n,j] == 3
                numFRight = numFlux.(obj, uQ[:,j],uQ[:,j+1]);
                numFLeft = numFlux.(obj, uQ[:,j-1],uQ[:,j]);
                uQNew[:,j] = uQ[:,j] -dt/dx*(numFRight-numFLeft);
                for i = 1:N
                    uNew[i,j] = IntegralVec(obj.q,uQNew[:,j].*obj.basis.PhiQuad[:,i]*0.5,-1.0,1.0);
                end
            elseif obj.refine[n,j] == 1 || obj.refine[n,j] == 2
                numFRight = numFlux.(obj, Ubb(obj.closureF,EvalAtQuad(obj.basisF,v[:,j])),Ubb(obj.closure,EvalAtQuad(obj.basisF,v[:,j+1])));
                numFLeft = numFlux.(obj, Ubb(obj.closureF,EvalAtQuad(obj.basisF,v[:,j-1])),Ubb(obj.closure,EvalAtQuad(obj.basisF,v[:,j])));
                for i = 1:N
                    uNew[i,j] = u[i,j]-dt/dx*IntegralVec(obj.qF,(numFRight-numFLeft).*obj.basisF.PhiQuad[:,i]*0.5,-1.0,1.0);
                end
            end
        end

        UpdateRefinement(obj,uNew,n+1);

        for j = 1:Nx
            if obj.refine[n+1,j] == 1 || obj.refine[n+1,j] == 2 || obj.refine[n+1,j] == 3
                if obj.refine[n,j] == 0 || obj.refine[n,j] == 3
                    vQ = -log.(obj.settings.uPlus - uQNew[:,j]) + log.(uQNew[:,j]-obj.settings.uMinus);
                    uQF = Ubb(obj.closure,obj.L*vQ);
                    uNew[:,j] = ComputeMomentsWithQuad(obj.closureF,uQF)
                end
                # (can be left out for collocation version) project to IPM solution
                vNew[:,j] = Solve(obj.closureF,uNew[:,j],v[:,j]);
                v[:,j] = vNew[:,j];
                u[:,j] = ComputeMoments(obj.closureF,vNew[:,j]);
            end
            if obj.refine[n+1,j] == 0 || obj.refine[n+1,j] == 2 || obj.refine[n+1,j] == 3
                if obj.refine[n,j] == 0 || obj.refine[n,j] == 3
                    uQ[:,j] = uQNew[:,j];
                else
                    vNew[:,j] = Solve(obj.closureF,uNew[:,j],v[:,j]);
                    v[:,j] = vNew[:,j];
                    uQ[:,j] = Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j]));
                end
                u[:,j] = uNew[:,j];
            end  
        end
        t = t+dt;
    end

    for j = 2:(Nx-1)
        uQ[:,j] = uQNew[:,j];
        if obj.refine[end-1,j] == 0 || obj.refine[end-1,j] == 3
            vQ = -log.(obj.settings.uPlus - uQ[:,j]) + log.(uQ[:,j]-obj.settings.uMinus);
            uQF = Ubb(obj.closure,obj.L*vQ);
            uNew[:,j] = ComputeMomentsWithQuad(obj.closureF,uQF)
        end
        vNew[:,j] = Solve(obj.closureF,uNew[:,j],v[:,j]);
        v[:,j] = vNew[:,j];
        #u[:,j] = ComputeMoments(obj.closureF,vNew[:,j]);
    end

    # return end time and solution
    return t, u, v;

end

function UpdateRefinement(obj::IPMSolverAdaptive,u::Array{Float64,2},n::Int64)
    for j = 1:obj.settings.Nx # find fine cells
        uL2 = IntegralVec(obj.q,EvalAtQuad(obj.basis,u[:,j]).^2,-1.0,1.0);
        indicator = u[end,j]^2*IntegralVec(obj.q,obj.basis.PhiQuad[:,end].^2,-1.0,1.0) / uL2;
        indicator += u[end-1,j]^2*IntegralVec(obj.q,obj.basis.PhiQuad[:,end-1].^2,-1.0,1.0) / uL2;
        if indicator > obj.settings.barrierRefine
            obj.refine[n,j] = 1;
        end
    end
    for j = 2:(obj.settings.Nx-1) # make all coarse cells neighboring fine cells type 2
        if obj.refine[n,j] == 1  
            if obj.refine[n,j-1] == 0
                obj.refine[n,j-1] = 2;
            end
            if obj.refine[n,j+1] == 0
                obj.refine[n,j+1] = 2;
            end
        end
    end
    for j = 2:(obj.settings.Nx-1) # make all coarse cells neighboring type 2 cells type 3
        if obj.refine[n,j] == 2  
            if obj.refine[n,j-1] == 0
                obj.refine[n,j-1] = 3;
            end
            if obj.refine[n,j+1] == 0
                obj.refine[n,j+1] = 3;
            end
        end
    end

end

function f(obj::IPMSolverAdaptive,u)
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