__precompile__
using PyPlot

struct Plotting
    settings::Settings;
    basis::Basis;
    xQuad::Array{Float64,1};
    x::Array{Float64,1}
    tEnd::Float64;
    q::Quadrature;
    euler::Euler;

    function Plotting(settings::Settings,basis::Basis,quadrature::Quadrature,tEnd::Float64)
        euler = Euler(settings);
        new(settings,basis,quadrature.xi,settings.x,tEnd,quadrature,euler);
    end
end

function PlotInX(obj::Plotting,s::Int,u::Array{Float64,3},xi::Array{Float64,1})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxFine = 1000;
    xFine = linspace(obj.settings.a,obj.settings.b,NxFine)
    uExact = zeros(NxFine);
    uPlot = zeros(Nx);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    for k = 1:length(xi)
        for j = 1:Nx
            tmp = Eval(obj.basis,u[:,s,j],xi[k]);
            uPlot[j] = tmp[1];
        end
        #for j = 1:NxFine
            #uExact[j] = obj.settings.solutionExact(obj.tEnd,xFine[j],xi[k])[1];
        #end
        ax[:plot](obj.x,uPlot, "r--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
        #ax[:plot](xFine,uExact, "k-", linewidth=2, label=L"$u_{ex}$", alpha=0.6)
    end
    #ylimMinus = -2.0;
    #ylimPlus = 15.0
    #ax[:set_ylim]([ylimMinus,ylimPlus])
    ax[:set_xlim]([obj.settings.a,obj.settings.b])
    ax[:set_xlabel]("x", fontsize=20);

end

function PlotInXi(obj::Plotting,s::Int,u::Array{Float64,3},index::Int)
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxiFine = 100;
    xiFine = linspace(-1,1,NxiFine)
    uExact = zeros(NxiFine);
    uPlot = zeros(NxiFine);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    if s == 1
        title(L"$\rho$")
    elseif s == 2
        title(L"$\rho u$")
    elseif s == 3
        title(L"$\rho v$")
    else
        title(L"$\rho e$")
    end
    for k = 1:NxiFine
        tmp = Eval(obj.basis,u[:,s,index],xiFine[k]);
        uPlot[k] = tmp[1];
        data_rho,data_u,data_P,data_e = obj.settings.solutionExact(obj.tEnd,obj.x[index],xiFine[k]);
        if s == 1
            uExact[k] = data_rho[1];
        elseif s == 2
            uExact[k] = data_rho[1]*data_u[1];
        else
            uExact[k] = data_rho[1]*data_e[1];
        end
    end

    ax[:plot](xiFine,uPlot, "r--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
    ax[:plot](xiFine,uExact, "k-", linewidth=2, label=L"$u_{ex}$", alpha=0.6)
    ax[:set_xlim]([-1.0,1.0])
    ax[:set_xlabel](L"\xi", fontsize=20);

end

function ComparePlotInXi(obj::Plotting,s::Int,u::Array{Float64,3},v::Array{Float64,3},index::Int)
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxiFine = 100;
    xiFine = linspace(-1,1,NxiFine)
    uExact = zeros(NxiFine);
    uPlot = zeros(NxiFine);
    vPlot = zeros(NxiFine);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    if s == 1
        title(L"$\rho$")
    elseif s == 2
        title(L"$\rho u$")
    elseif s == 3
        title(L"$\rho v$")
    else
        title(L"$\rho e$")
    end
    for k = 1:NxiFine
        tmp = Eval(obj.basis,u[:,s,index],xiFine[k]);
        uPlot[k] = tmp[1];
        tmp = Eval(obj.basis,v[:,s,index],xiFine[k]);
        vPlot[k] = tmp[1];
        data_rho,data_u,data_P,data_e = obj.settings.solutionExact(obj.tEnd,obj.x[index],xiFine[k]);
        if s == 1
            uExact[k] = data_rho[1];
        elseif s == 2
            uExact[k] = data_rho[1]*data_u[1];
        else
            uExact[k] = data_rho[1]*data_e[1];
        end
    end

    ax[:plot](xiFine,uPlot, "r--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
    ax[:plot](xiFine,vPlot, "b-.", linewidth=2, label=L"$u_{N}$", alpha=0.6)
    ax[:plot](xiFine,uExact, "k-", linewidth=2, label=L"$u_{ex}$", alpha=0.6)
    ax[:set_xlim]([-1.0,1.0])
    ax[:set_xlabel](L"\xi", fontsize=15);

end

function PlotExpectedValue(obj::Plotting,s::Int,u::Array{Float64,4})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    Ny = obj.settings.Ny;
    NxFine = 1000;
    xFine = linspace(obj.settings.a,obj.settings.b,NxFine)
    x = linspace(obj.settings.a+0.5*obj.settings.dx,obj.settings.b-0.5*obj.settings.dx,Nx)
    uExact = zeros(NxFine);
    varExact = zeros(NxFine);
    uPlot = zeros(Nx);
    varPlot = zeros(Nx);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    if s == 1
        title(L"$\rho$", fontsize=20)
    elseif s == 2
        title(L"$\rho u$", fontsize=20)
    elseif s == 3
        title(L"$\rho v$", fontsize=20)
    else
        title(L"$\rho e$", fontsize=20)
    end
    for j = 1:Nx
        i = Int64(round(0.5*obj.settings.Ny));
        uPlot[j] = u[1,s,j,i];
        varPlot[j] = u[2:end,s,j,i]'u[2:end,s,j,i];
    end
    varMax = maximum(varPlot);
    expMax = maximum(uPlot);
    qFine = Quadrature(200,"Gauss")
    exactState = zeros(NxFine,qFine.Nq);
    for k = 1:qFine.Nq
        #data_rho,data_u,data_P,data_e = analytic_sod(obj.tEnd,obj.settings.x0,obj.settings.rhoL,obj.settings.pL,obj.settings.uL,obj.settings.rhoR,
        #                                obj.settings.pR,obj.settings.uR,obj.settings.gamma,xFine);
        data_rho,data_u,data_P,data_e = obj.settings.solutionExact(obj.tEnd,xFine,qFine.xi[k]);
        if s == 1
            exactState[:,k] = data_rho;
        elseif s == 2
            exactState[:,k] = data_rho.*data_u;
        else
            exactState[:,k] = data_rho.*data_e;
        end
    end
    for j = 1:NxFine
        uExact[j] = IntegralVec(qFine, exactState[j,:]*0.5,-1.0,1.0);
        varExact[j] = IntegralVec(qFine, (exactState[j,:]-uExact[j]).^2*0.5,-1.0,1.0)
    end
    ax[:plot](x,uPlot, "r-", linewidth=2, label=L"$E[u_{N}]$", alpha=0.6)
    #ax[:plot](xFine,uExact, "r-", linewidth=2, label=L"$E[u_{ex}]$", alpha=0.6)
    ylabel("Expectation Value", fontsize=20,color="red")
    ax2 = ax[:twinx]() # Create another axis on top of the current axis
    #font2 = ["color"=>"blue"]
    ylabel("Variance", fontsize=20,color="blue")
    ax2[:plot](x,varPlot, "b-", linewidth=2, label=L"$Var[u_{N}]$", alpha=0.6)
    setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
    setp(ax[:get_yticklabels](),color="red")
    setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
    setp(ax[:get_yticklabels](),color="red")
    ax[:set_xlim]([obj.settings.a,obj.settings.b])
    ax[:set_xlabel]("x", fontsize=20);  
    ax[:tick_params]("both",labelsize=20) 
    ax2[:tick_params]("both",labelsize=20)
    fig[:canvas][:draw]() # Update the figure
    savefig("results/XCutState$(s)Filter$(obj.settings.filterType)Nx$(Nx)N$(obj.settings.N)tEnd$(obj.settings.tEnd)Sigma$(obj.settings.sigma).png")
end

function PlotExpectedValueYSlice(obj::Plotting,s::Int,u::Array{Float64,4})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    Ny = obj.settings.Ny;
    NxFine = 1000;
    xFine = linspace(obj.settings.c,obj.settings.d,NxFine)
    x = linspace(obj.settings.c+0.5*obj.settings.dy,obj.settings.d-0.5*obj.settings.dy,Ny)
    uExact = zeros(NxFine);
    varExact = zeros(NxFine);
    uPlot = zeros(Ny);
    varPlot = zeros(Ny);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    if s == 1
        title(L"$\rho$", fontsize=20)
    elseif s == 2
        title(L"$\rho u$", fontsize=20)
    elseif s == 3
        title(L"$\rho v$", fontsize=20)
    else
        title(L"$\rho e$", fontsize=20)
    end
    for j = 1:Ny
        i = Int64(round(0.5*obj.settings.Nx));
        uPlot[j] = u[1,s,i,j];
        varPlot[j] = u[2:end,s,i,j]'u[2:end,s,i,j];
    end
    varMax = maximum(varPlot);
    expMax = maximum(uPlot);
    qFine = Quadrature(200,"Gauss")
    exactState = zeros(NxFine,qFine.Nq);
    for k = 1:qFine.Nq
        #data_rho,data_u,data_P,data_e = analytic_sod(obj.tEnd,obj.settings.x0,obj.settings.rhoL,obj.settings.pL,obj.settings.uL,obj.settings.rhoR,
        #                                obj.settings.pR,obj.settings.uR,obj.settings.gamma,xFine);
        data_rho,data_u,data_P,data_e = obj.settings.solutionExact(obj.tEnd,xFine,qFine.xi[k]);
        if s == 1
            exactState[:,k] = data_rho;
        elseif s == 2
            exactState[:,k] = data_rho.*data_u;
        else
            exactState[:,k] = data_rho.*data_e;
        end
    end
    for j = 1:NxFine
        uExact[j] = IntegralVec(qFine, exactState[j,:]*0.5,-1.0,1.0);
        varExact[j] = IntegralVec(qFine, (exactState[j,:]-uExact[j]).^2*0.5,-1.0,1.0)
    end
    ax[:plot](x,uPlot, "r-", linewidth=2, label=L"$E[u_{N}]$", alpha=0.6)
    #ax[:plot](xFine,uExact, "r-", linewidth=2, label=L"$E[u_{ex}]$", alpha=0.6)
    ylabel("Expectation Value", fontsize=20,color="red")
    ax2 = ax[:twinx]() # Create another axis on top of the current axis
    #font2 = ["color"=>"blue"]
    ylabel("Variance",color="blue")
    ax2[:plot](x,varPlot, "b-", linewidth=2, label=L"$Var[u_{N}]$", alpha=0.6)
    setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
    setp(ax[:get_yticklabels](),color="red")
    #ax2[:plot](xFine,varExact, "b-", linewidth=2, label=L"$V[u_{ex}]$", alpha=0.6)
    #ylimMinus = -0.5;
    #ylimPlus = 16.0
    #ax[:set_ylim]([uExact[end]-0.1,uExact[1]+0.1])
    #ax[:set_ylim]([uExact[end]-0.1,uExact[1]+0.1])
    setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
    setp(ax[:get_yticklabels](),color="red")
    ax[:set_xlim]([obj.settings.c,obj.settings.d])
    ax[:set_xlabel]("x", fontsize=20);  
    ax[:tick_params]("both",labelsize=20) 
    ax2[:tick_params]("both",labelsize=20)
    fig[:canvas][:draw]() # Update the figure
    savefig("results/YCutState$(s)Filter$(obj.settings.filterType)Nx$(Nx)N$(obj.settings.N)tEnd$(obj.settings.tEnd)Sigma$(obj.settings.sigma).png")
end

function Plot2DExpectedValue(obj::Plotting,s::Int,u::Array{Float64,4})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    Ny = obj.settings.Ny;
    NxFine = 1000;
    xFine = linspace(obj.settings.c,obj.settings.d,NxFine)
    x = linspace(obj.settings.c+0.5*obj.settings.dy,obj.settings.d-0.5*obj.settings.dy,Ny)
    uExact = zeros(NxFine);
    varExact = zeros(NxFine);
    uPlot = zeros(Ny);
    varPlot = zeros(Ny);

    # start plot
    fig, ax = subplots(figsize=(10.5, 8), dpi=100)
    if s == 1
        title(L"$\rho$", fontsize=20)
    elseif s == 2
        title(L"$\rho u$", fontsize=20)
    elseif s == 3
        title(L"$\rho v$", fontsize=20)
    else
        title(L"$\rho e$", fontsize=20)
    end
    
    pcolormesh(u[1,s,:,:]')
    colorbar()
    savefig("results/ExpState$(s)Filter$(obj.settings.filterType)Nx$(Nx)N$(obj.settings.N)tEnd$(obj.settings.tEnd)Sigma$(obj.settings.sigma).png")
end