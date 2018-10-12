__precompile__
include("quadrature.jl")
include("Basis.jl")
include("euler.jl")
include("mesh.jl")
include("closure.jl")

using Plots

using ProgressMeter

struct Solver
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # quadrature
    q::Quadrature;

    # Solver settings
    settings::Settings;

    # spatial basis functions
    basis::Basis;

    # filter strength
    lambdaFilter::Float64;

    # problem type
    euler::Euler;

    # Euler 2D closure
    closure::Closure;

    # mesh
    mesh::Mesh;

    # mesh with holes
    flowVolume::Array{Float64,2};

    # preallocated memory for fluxes
    fluxNorth::Array{Float64,2};
    fluxEast::Array{Float64,2};
    fluxSouth::Array{Float64,2};
    fluxWest::Array{Float64,2};

    # constructor
    function Solver(settings)
        x = linspace(settings.a,settings.b,settings.Nx)
        q = Quadrature(settings.Nq,"Gauss");
        basis = Basis(q,settings);
        euler = Euler(settings);
        mesh = Mesh(settings);
        closure = Closure(settings,basis,q)

        flowVolume = zeros(settings.Nx,settings.Ny);

        xM = zeros(2); xM[1] = 0.5; xM[2] = 0.35;
        flowVolume = InsertRectangularObstacle(mesh,xM,0.2,0.35,flowVolume);

        flowVolume = AddWallBoundary(mesh,"west", flowVolume);
        flowVolume = AddWallBoundary(mesh,"east", flowVolume);
        flowVolume = AddWallBoundary(mesh,"north", flowVolume);

        fluxNorth = zeros(settings.Nq,settings.states);
        fluxEast = zeros(settings.Nq,settings.states);
        fluxSouth = zeros(settings.Nq,settings.states);
        fluxWest = zeros(settings.Nq,settings.states);

        new(x,q,settings,basis,s.lambdaFilter,euler,closure,mesh,flowVolume,fluxNorth,fluxEast,fluxSouth,fluxWest);
    end
end

function numFluxEuler(obj::Solver,u,v,nUnit::Array{Float64,1},n::Array{Float64,1})
    # numerical HLL Flux
    uU,vU,aU = SpeedFromCons(obj.euler,u);
    uV,vV,aV = SpeedFromCons(obj.euler,v);

    uUProjected = nUnit[1]*uU + nUnit[2]*vU;
    uVProjected = nUnit[1]*uV + nUnit[2]*vV;
    
    lambdaMin = uUProjected-aU;
    lambdaMax = uVProjected+aV;

    if lambdaMin >= 0
        y = physicalFlux(obj.euler,u)*n;
    elseif lambdaMax <= 0
        y = physicalFlux(obj.euler,v)*n;
    else
        factor = 1.0/(lambdaMax-lambdaMin);
        flux1 = physicalFlux(obj.euler,u);
        flux2 = physicalFlux(obj.euler,v);
        y = factor*(lambdaMax*flux1*n-lambdaMin*flux2*n+lambdaMax*lambdaMin*(v-u)*norm(n));
    end
    return y;
end

function rhs(obj::Solver,u::Array{Float64,2},uN::Array{Float64,2},uE::Array{Float64,2},uS::Array{Float64,2},uW::Array{Float64,2},i::Int,j::Int)

    faceNormals = obj.mesh.faceNormals[i,j,:,:];
    for k = 1:obj.settings.Nq
        if obj.flowVolume[i-1,j] == 1 # n = [1 0]
            uG = u[k,:];
            uG[2] = -uG[2];
            obj.fluxWest[k,:] = numFluxEuler(obj,uG,u[k,:],obj.mesh.faceUnitNormals[i,j,1,:],faceNormals[1,:]);
        else
            obj.fluxWest[k,:] = numFluxEuler(obj,uW[k,:],u[k,:],obj.mesh.faceUnitNormals[i,j,1,:],faceNormals[1,:]);
        end

        if obj.flowVolume[i+1,j] == 1 # n = [-1 0]
            uG = u[k,:];
            uG[2] = -uG[2];
            obj.fluxEast[k,:] = numFluxEuler(obj,u[k,:],uG,obj.mesh.faceUnitNormals[i,j,3,:],faceNormals[3,:]);
        else
            obj.fluxEast[k,:] = numFluxEuler(obj,u[k,:],uE[k,:],obj.mesh.faceUnitNormals[i,j,3,:],faceNormals[3,:]);
        end

        if obj.flowVolume[i,j+1] == 1 # n = [0 -1]
            uG = u[k,:];
            uG[3] = -uG[3];
            obj.fluxNorth[k,:] = numFluxEuler(obj,u[k,:],uG,obj.mesh.faceUnitNormals[i,j,2,:],faceNormals[2,:]);
        else
            obj.fluxNorth[k,:] = numFluxEuler(obj,u[k,:],uN[k,:],obj.mesh.faceUnitNormals[i,j,2,:],faceNormals[2,:]);
        end
        if obj.flowVolume[i,j-1] == 1 # n = [0 1]
            uG = u[k,:];
            uG[3] = -uG[3];
            obj.fluxSouth[k,:] = numFluxEuler(obj,uG,u[k,:],obj.mesh.faceUnitNormals[i,j,4,:],faceNormals[4,:]);
        else
            obj.fluxSouth[k,:] = numFluxEuler(obj,uS[k,:],u[k,:],obj.mesh.faceUnitNormals[i,j,4,:],faceNormals[4,:]);
        end
    end

    totalFlux = (obj.fluxNorth-obj.fluxSouth)+(obj.fluxEast-obj.fluxWest);

    return (1.0/obj.mesh.cellVolumes[i,j])*ComputeMoments(obj.basis,totalFlux*0.5);
end

function SetupIC(obj::Solver)
    u = zeros(obj.settings.N,obj.settings.states,obj.settings.Nx,obj.settings.Ny);
    for j = 1:obj.settings.Nx
        for l = 1:obj.settings.Ny
            for s = 1:obj.settings.states
                for i = 1:obj.settings.N
                    # check if obstacle
                    if obj.flowVolume[j,l] == 1
                        u[i,s,j,l] = NaN;
                    else
                        u[i,s,j,l] = IntegralVec(obj.q, obj.settings.IC(s,obj.mesh.cellMidPoints[j,l,:],obj.q.xi).*obj.basis.PhiQuad[:,i]*0.5);
                    end
                end
            end
        end
    end
    return u;
end

function Solve(obj::Solver)
    t = 0.0;
    dt = obj.settings.dt;
    dx = obj.x[2]-obj.x[1];
    tEnd = obj.settings.tEnd;
    Nx = obj.settings.Nx;
    Ny = obj.settings.Ny;
    N = obj.settings.N; # number of moments

    var = zeros(Nx,Ny);

    # Set up initial condition
    println("Setting up IC");
    u = SetupIC(obj);

    v = zeros(size(u));
    vNew = zeros(size(u));
    for j = 1:Ny
        for i = 1:Nx
            if obj.flowVolume[i,j] != 1
                v[1,4,i,j] = -1.0; v[1,1,i,j] = 1.0; # ensure log(-v4) is real
                vNew[:,:,i,j] = Solve(obj.closure,u[:,:,i,j],v[:,:,i,j]);
                v[:,:,i,j] = vNew[:,:,i,j];
                u[:,:,i,j] = ComputeMoments(obj.closure,vNew[:,:,i,j]);
            end
        end
    end
    uNew = deepcopy(u);

    # check CFL
    cfl = getCFL(obj.euler,u[1,:,1,1],dt,dx);
    println("CFL is ",cfl);    

    Nt = round(tEnd/dt);
    Nt::Int = round(tEnd/dt);
    prog = Progress(Nt,1)

    #pyplot(size=(800,600))

    clibraries()
    clibrary(:colorcet)

    # choose state for plotting. 1 = density etc
    plotState = 1;

    #create initial gif
    for j = 1:Ny
        for i = 1:Nx
            var[i,j] = u[2:end,plotState,i,j]'u[2:end,plotState,i,j];
        end
    end
    #q = Plots.heatmap(u[1,plotState,:,:],colorbar=true,color=:viridis,ratio=:equal)
    p = Plots.heatmap(var',colorbar=true,clim=(0.0,0.06),color=:viridis,ratio=:equal)

    # time loop    
    #@gif 
    for n = 1:Nt
    #for n = 1:Nt

        # Update time by dt
        for j = 2:(Ny-1)
            for i = 2:(Nx-1)

                if obj.flowVolume[i,j] == 1
                    continue;
                end

                rhsVals = rhs(obj, UKin(obj.closure,EvalAtQuad(obj.basis,v[:,:,i,j])), UKin(obj.closure,EvalAtQuad(obj.basis,v[:,:,i,j+1])),UKin(obj.closure,EvalAtQuad(obj.basis,v[:,:,i+1,j])),
                    UKin(obj.closure,EvalAtQuad(obj.basis,v[:,:,i,j-1])),UKin(obj.closure,EvalAtQuad(obj.basis,v[:,:,i-1,j])),i,j);
                uNew[:,:,i,j] = u[:,:,i,j]-dt*rhsVals;
            end
        end

        for j = 2:(Ny-1)
            for i = 2:(Nx-1)
                if obj.flowVolume[i,j] == 1
                    continue;
                end
                vNew[:,:,i,j] = Solve(obj.closure,uNew[:,:,i,j],v[:,:,i,j]);
                v[:,:,i,j] = vNew[:,:,i,j];
                u[:,:,i,j] = ComputeMoments(obj.closure,vNew[:,:,i,j]);
            end
        end
        t = t+dt;
        # expectation of plotState
        
        # variance of plotState
        #for j = 1:Ny
        #    for i = 1:Nx
        #        var[i,j] = u[2:end,plotState,i,j]'u[2:end,plotState,i,j];
        #    end
        #end
        #q = Plots.heatmap!(q[1],u[1,plotState,:,:],colorbar=true,color=:viridis,ratio=:equal) #color=:coolwarm,:rainbow, :viridis ,clim=(0.0,0.04)
        #p = Plots.heatmap!(p[1],var',colorbar=true,clim=(0.0,0.06),color=:viridis,ratio=:equal) #color=:coolwarm,:rainbow, :viridis ,clim=(0.0,0.04)

        next!(prog) # update progress bar
    end #every 2

    println("tEnd is ",t)

    # return end time and solution
    return t, u;

end
