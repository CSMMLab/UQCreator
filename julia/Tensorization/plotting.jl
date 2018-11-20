__precompile__
using PyPlot

struct Plotting
    settings::Settings;
    basis::Basis;
    xQuad::Array{Float64,1};
    x::Array{Float64,1}
    tEnd::Float64;
    q::Quadrature;
    closure::Closure;

    function Plotting(settings::Settings,basis::Basis,quadrature::Quadrature,closure::Closure,tEnd::Float64)
        new(settings,basis,quadrature.xi,settings.x,tEnd,quadrature,closure);
    end
end

function PlotInX(obj::Plotting,v::Array{Float64,2},xi::Array{Float64,1})
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
            tmp = Ubb(obj.closure,Eval(obj.basis,v[:,j],xi[k]));
            uPlot[j] = tmp[1];
        end
        for j = 1:NxFine
            uExact[j] = obj.settings.solutionExact(obj.tEnd,xFine[j],xi[k])[1];
        end
        ax[:plot](obj.x,uPlot, "r--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
        ax[:plot](xFine,uExact, "k-", linewidth=2, label=L"$u_{ex}$", alpha=0.6)
    end
    ylimMinus = -2.0;
    ylimPlus = 15.0
    ax[:set_ylim]([ylimMinus,ylimPlus])
    ax[:set_xlim]([obj.settings.a,obj.settings.b])
    ax[:set_xlabel]("x", fontsize=20);

end

function PlotInXi(obj::Plotting,v::Array{Float64,2},index::Int)
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxiFine = 300;
    xiFine = linspace(-1,1,NxiFine)
    uExact = zeros(NxiFine);
    uPlot = zeros(NxiFine);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    for k = 1:NxiFine
        tmp = Ubb(obj.closure,Eval(obj.basis,v[:,index],xiFine[k]));
        uPlot[k] = tmp[1];
        uExact[k] = obj.settings.solutionExact(obj.tEnd,obj.x[index],xiFine[k])[1];
    end

    ax[:plot](xiFine,uPlot, "r--", linewidth=2, label=L"$u_{N}$", alpha=0.6)
    ax[:plot](xiFine,uExact, "k-", linewidth=2, label=L"$u_{ex}$", alpha=0.6)
    ylimMinus = -2.0;
    ylimPlus = 15.0
    ax[:set_ylim]([ylimMinus,ylimPlus])
    ax[:set_xlim]([-1.0,1.0])
    ax[:set_xlabel](L"\xi", fontsize=20);

end

function PlotExpectedValue(obj::Plotting,v::Array{Float64,2})
    Nq = obj.settings.Nq;
    Nx = obj.settings.Nx;
    NxFine = 1000;
    xFine = linspace(obj.settings.a,obj.settings.b,NxFine)
    uExact = zeros(NxFine);
    varExact = zeros(NxFine);
    uPlot = zeros(Nx);
    varPlot = zeros(Nx);

    # start plot
    fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    for j = 1:Nx
        uVals = Ubb(obj.closure,EvalAtQuad(obj.basis,v[:,j]));
        uPlot[j] = IntegralVec(obj.q,uVals*0.5,-1.0,1.0);
        varPlot[j] = IntegralVec(obj.q,0.5*(uVals-uPlot[j]).^2,-1.0,1.0);

    end
    qFine = Quadrature(200,"Gauss");
    for j = 1:NxFine
        uExact[j] = Integral(qFine, xi-> obj.settings.solutionExact(obj.tEnd,xFine[j],xi)*0.5,-1.0,1.0);
        varExact[j] = Integral(qFine, xi-> (obj.settings.solutionExact(obj.tEnd,xFine[j],xi)-uExact[j]).^2*0.5,-1.0,1.0);
    end
    ax[:plot](obj.x,uPlot, "r--", linewidth=2, label=L"$E[u_{N}]$", alpha=0.6)
    ax[:plot](xFine,uExact, "k-", linewidth=2, label=L"$E[u_{ex}]$", alpha=0.6)
    ax[:plot](obj.x,0.5*varPlot, "r:", linewidth=2, label=L"$V[u_{N}]$", alpha=0.6)
    ax[:plot](xFine,0.5*varExact, "k:", linewidth=2, label=L"$V[u_{ex}]$", alpha=0.6)
    ylimMinus = -0.5;
    ylimPlus = 16.0
    ax[:set_ylim]([ylimMinus,ylimPlus])
    ax[:set_xlim]([obj.settings.a,obj.settings.b])
    ax[:set_xlabel]("x", fontsize=20);
end

function L2Error(obj::Plotting,v::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    error = 0;
    qFine = Quadrature(300,"Gauss");
    basisFine = Basis(qFine,obj.settings);
    for j = 1:Nx
        uVals = Ubb(obj.closure,EvalAtQuad(basisFine,v[:,j]));
        uExVals = obj.settings.solutionExact(t,x[j],qFine.xi);
        error = error + obj.settings.dx*IntegralVec(qFine, 0.5*( uExVals-uVals).^2,-1.0,1.0)
    end
    return sqrt(error);
end

function L1Error(obj::Plotting,v::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    error = 0.0;
    qFine = Quadrature(300,"Gauss");
    basisFine = Basis(qFine,obj.settings);
    for j = 1:Nx
        uVals = Ubb(obj.closure,EvalAtQuad(basisFine,v[:,j]));
        uExVals = obj.settings.solutionExact(t,x[j],qFine.xi);
        error = error + obj.settings.dx*IntegralVec(qFine, 0.5*abs.( uExVals-uVals),-1.0,1.0)
        if x[j] < 1.0 && j == 1
            #println("cell ",j,": error is ",error)
        end
        if j == 1
            #println(uExVals)
            #fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
            #ax[:plot](obj.q.xi,uVals, "r--", linewidth=2, label=L"IPM", alpha=0.6)
            #ax[:plot](obj.q.xi,uExVals, "k-", linewidth=2, label=L"exact", alpha=0.6)
        end
    end
    return error;
end

function L2ErrorExpVar(obj::Plotting,v::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    errorExp = 0.0;
    errorVar = 0.0;
    qFine = Quadrature(300,"Gauss");
    basisFine = Basis(qFine,obj.settings);
    for j = 1:Nx
        uExVals = obj.settings.solutionExact(t,x[j],qFine.xi);
        uVals = Ubb(obj.closure,EvalAtQuad(basisFine,v[:,j]));
        expN = IntegralVec(qFine,uVals*0.5,-1.0,1.0);
        varN = IntegralVec(qFine,0.5*(uVals-expN).^2,-1.0,1.0);
        expEx = IntegralVec(qFine, uExVals*0.5,-1.0,1.0)
        varEx = IntegralVec(qFine, (uExVals-expEx).^2*0.5,-1.0,1.0);
        errorExp = errorExp + obj.settings.dx*(expEx-expN)^2;
        errorVar = errorVar + obj.settings.dx*(varEx-varN)^2;
    end
    return sqrt(errorExp),sqrt(errorVar);
end

function L1ErrorExpVar(obj::Plotting,v::Array{Float64,2},t::Float64)
    Nx = obj.settings.Nx;
    x = obj.settings.x;
    errorExp = 0.0;
    errorVar = 0.0;
    qFine = Quadrature(300,"Gauss");
    basisFine = Basis(qFine,obj.settings);
    for j = 1:Nx
        uExVals = obj.settings.solutionExact(t,x[j],qFine.xi);
        uVals = Ubb(obj.closure,EvalAtQuad(basisFine,v[:,j]));
        expN = IntegralVec(qFine,uVals*0.5,-1.0,1.0);
        varN = IntegralVec(qFine,0.5*(uVals-expN).^2,-1.0,1.0);
        expEx = IntegralVec(qFine, uExVals*0.5,-1.0,1.0)
        varEx = IntegralVec(qFine, (uExVals-expEx).^2*0.5,-1.0,1.0);
        errorExp = errorExp + obj.settings.dx*abs(expEx-expN);
        errorVar = errorVar + obj.settings.dx*abs(varEx-varN);
    end
    return errorExp,errorVar;
end