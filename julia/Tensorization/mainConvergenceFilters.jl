using Revise
include("settings.jl")
include("IPMSolver.jl")
include("plotting.jl")

function L1Error(s::Settings,N::Int)
    s.N = N;
    s.Nq = max(5*N,40);
    solver = IPMSolver(s);
    #@time tEnd, u, v = Solve(solver);
    runTime = @elapsed tEnd, u, v = Solve(solver);

    plotSolution = Plotting(s,solver.basis,solver.q,solver.closure,tEnd);
    l2err = L2Error(plotSolution,v,tEnd);
    l1err = L1Error(plotSolution,v,tEnd);
    errorExp,errorVar = L1ErrorExpVar(plotSolution,v,tEnd);
    return l2err,l1err,errorExp,errorVar,runTime;
end

s = Settings();
NMax = 45; # 45
NMin = 5; # 5
N = NMin:5:NMax;
errors = zeros(size(N));
errorsL1 = zeros(size(N));
errorsExp = zeros(size(N));
errorsVar = zeros(size(N));
runtimes = zeros(size(N));

s.filterType = "L1FilterClosure";
for i = 1:length(errors)
    (errors[i],errorsL1[i],errorsExp[i],errorsVar[i],runtimes[i]) = L1Error(s,N[i]);
    println("Error IPM for N = ",N[i]-1, " is " , errorsL1[i]);
end

fig, ax = subplots()
ax[:plot](N,errors, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,(errors[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("L2 error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsL1, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,(errorsL1[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("L1 error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsExp, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,(errorsExp[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("Exp error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsVar, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,(errorsVar[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("Var error", fontsize=20);

writedlm("NValsIPM.txt", N)

writedlm("errorsIPMVarNx$(s.Nx).txt", errorsVar)

writedlm("errorsIPMExpNx$(s.Nx).txt", errorsExp)

writedlm("errorsIPML1Nx$(s.Nx).txt", errorsL1)

writedlm("errorsIPMNx$(s.Nx).txt", errors)

writedlm("runtimesIPMNx$(s.Nx).txt", runtimes)

println("test finished")