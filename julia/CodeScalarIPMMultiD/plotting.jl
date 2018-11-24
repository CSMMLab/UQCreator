__precompile__
using PyPlot
using QuantEcon: meshgrid



struct Plotting
    settings::Settings;
    basis::Basis;
    xQuad::Array{Float64,1};
    x::Array{Float64,1}
    tEnd::Float64;
    q::Quadrature;
    closure;

    function Plotting(settings::Settings,basis::Basis,quadrature::Quadrature,closure,tEnd::Float64)
        new(settings,basis,quadrature.xi,settings.x,tEnd,quadrature,closure);
    end
end

function PlotInXiV(obj::Plotting,v::Array{Float64,2},index::Int)
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxiFine = 40; # 300
    xiFine = linspace(-1,1,NxiFine)
    uPlot = zeros(NxiFine,NxiFine);
    uEx = zeros(NxiFine,NxiFine);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    for k = 1:NxiFine
        for q = 1:NxiFine
            tmp = Ubb(obj.closure,Eval(obj.basis,v[:,index],xiFine[k],xiFine[q]));
            uPlot[k,q] = tmp[1];
            uEx[k,q] = obj.settings.solutionExact(obj.tEnd,obj.x[index],xiFine[k],xiFine[q])[1]
        end
    end

    xgrid, ygrid = meshgrid(xiFine, xiFine)
    surf(xgrid, ygrid, uPlot', cmap=ColorMap("jet"), alpha=0.7)

    fig = figure("solution IPM at x = $(obj.x[index])",figsize=(10,10))
    ax = fig[:add_subplot](2,1,1, projection = "3d")
    ax[:plot_surface](xgrid, ygrid, uPlot', rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    xlabel(L"$\xi_{0}$")
    ylabel(L"$\xi_{1}$")
    title("Solution $(obj.settings.closureType) IPM")

    subplot(212)
    ax = fig[:add_subplot](2,1,2)
    cp = ax[:contour](xgrid, ygrid, uPlot', colors="black", linewidth=2.0)
    ax[:clabel](cp, inline=1, fontsize=10)
    xlabel(L"$\xi_{0}$")
    ylabel(L"$\xi_{1}$")
    title("Solution Contours")
    tight_layout()

    fig = figure("solution exact at x = $(obj.x[index])",figsize=(10,10))
    ax = fig[:add_subplot](2,1,1, projection = "3d")
    ax[:plot_surface](xgrid, ygrid, uEx, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
    xlabel(L"$\xi_{0}$")
    ylabel(L"$\xi_{1}$")
    title("Exact Solution")

    subplot(212)
    ax = fig[:add_subplot](2,1,2)
    cp = ax[:contour](xgrid, ygrid, uEx, colors="black", linewidth=2.0)
    ax[:clabel](cp, inline=1, fontsize=10)
    xlabel(L"$\xi_{0}$")
    ylabel(L"$\xi_{1}$")
    title("Solution Contours")
    tight_layout()
end

function PlotExpectedValue(obj::Plotting,u::Array{Float64,2})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxFine = 500;
    xFine = linspace(obj.settings.a,obj.settings.b,NxFine)
    uExact = zeros(NxFine);
    varExact = zeros(NxFine);
    uPlot = zeros(Nx);
    varPlot = zeros(Nx);
    qFine = Quadrature(50,"Gauss")

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    for j = 1:Nx
        uPlot[j] = u[1,j];
        varPlot[j] = u[2:end,j]'u[2:end,j];
    end
    for j = 1:NxFine
        uExact[j] = Integral(qFine, (xi,eta)-> obj.settings.solutionExact(obj.tEnd,xFine[j],xi,eta)*0.25);
        varExact[j] = Integral(qFine, (xi,eta)-> 0.25*(obj.settings.solutionExact(obj.tEnd,xFine[j],xi,eta)-uExact[j]).^2);
    end

    ax[:plot](obj.x,uPlot, "k--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
    ylabel("Expectation Value", fontsize=20,color="red")
    ax[:plot](xFine,uExact, "r-", linewidth=2, alpha=0.6)
    ax2 = ax[:twinx]() # Create another axis on top of the current axis
    ylabel("Variance", fontsize=20,color="blue")
    ax2[:plot](obj.x,varPlot, "k:", linewidth=2, label="SG", alpha=0.6)
    #ax2[:set_position](new_position) # Position Method 2
    setp(ax2[:get_yticklabels](),color="blue") # Y Axis font formatting
    setp(ax[:get_yticklabels](),color="red")
    ax2[:plot](xFine,varExact, "b-", linewidth=2, alpha=0.6)
    #ylimMinus = -0.5;
    #ylimPlus = 16.0
    #ax[:set_ylim]([ylimMinus,ylimPlus])
    ax[:set_xlim]([obj.settings.a,obj.settings.b])
    ax[:set_xlabel]("x", fontsize=20);
    ax[:legend](loc="upper right", fontsize=20)
    ax[:tick_params]("both",labelsize=20) 
    ax2[:tick_params]("both",labelsize=20)
    fig[:canvas][:draw]() # Update the figure
    savefig("results/BurgersIC3ExpVarNx$(Nx)N$(obj.settings.N)tEnd$(obj.settings.tEnd)Sigma$(obj.settings.sigma0).png")
end

function L2Error(obj::Plotting,u::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    error = 0;
    for j = 1:Nx
        error = error + obj.settings.dx*Integral(obj.q, xi-> 0.25*( obj.settings.solutionExact(t,x[j],xi)-Eval(obj.basis,u[:,j],xi)).^2,-1.0,1.0)
    end
    return sqrt(error);
end

function L1Error(obj::Plotting,u::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    error = 0;
    for j = 1:Nx
        error = error + obj.settings.dx*Integral(obj.q, xi-> 0.25*abs.( obj.settings.solutionExact(t,x[j],xi)-Eval(obj.basis,u[:,j],xi)),-1.0,1.0)
    end
    return error;
end

function L2ErrorExpVar(obj::Plotting,u::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    errorExp = 0.0;
    errorVar = 0.0;
    for j = 1:Nx
        expN = u[1,j];
        varN = u[2:end,j]'u[2:end,j];
        expEx = Integral(obj.q, xi-> obj.settings.solutionExact(obj.tEnd,x[j],xi)*0.25,-1.0,1.0)
        varEx = Integral(obj.q, xi-> (obj.settings.solutionExact(obj.tEnd,x[j],xi)-expEx).^2*0.25,-1.0,1.0);
        errorExp = errorExp + obj.settings.dx*(expEx-expN)^2;
        errorVar = errorVar + obj.settings.dx*(varEx-varN)^2;
    end
    return sqrt(errorExp),sqrt(errorVar);
end

function L1ErrorExpVar(obj::Plotting,u::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    errorExp = 0.0;
    errorVar = 0.0;
    for j = 1:Nx
        expN = u[1,j];
        varN = u[2:end,j]'u[2:end,j];
        expEx = Integral(obj.q, xi-> obj.settings.solutionExact(obj.tEnd,x[j],xi)*0.25,-1.0,1.0)
        varEx = Integral(obj.q, xi-> (obj.settings.solutionExact(obj.tEnd,x[j],xi)-expEx).^2*0.25,-1.0,1.0);
        errorExp = errorExp + obj.settings.dx*abs(expEx-expN);
        errorVar = errorVar + obj.settings.dx*abs(varEx-varN);
    end
    return errorExp,errorVar;
end