using Revise
using TimeIt
include("settings.jl")
include("Solver.jl")
include("plotting.jl")

function L1Error(s::Settings,N::Int)
    s.N = N;
    s.Nq = 2*(N+5);
    solver = Solver(s);
    #runTime = @timeit tEnd, u = Solve(solver);
    runTime = @elapsed tEnd, u = Solve(solver);
    
    plotSolution = Plotting(s,solver.basis,solver.q,tEnd);
    l2err = L2Error(plotSolution,u,tEnd);
    l1err = L1Error(plotSolution,u,tEnd);
    errorExp,errorVar = L1ErrorExpVar(plotSolution,u,tEnd);
    return l2err,l1err,errorExp,errorVar,runTime;
end

s = Settings();
NMax = 70; # 70
NMin = 5; # 5
N = NMin:5:NMax;
errors = zeros(size(N));
errorsL1 = zeros(size(N));
errorsExp = zeros(size(N));
errorsVar = zeros(size(N));
runtimesL1 = zeros(size(N));

s.filterType = "L1FilterClosure";
for i = 1:length(errors)
    (errors[i],errorsL1[i],errorsExp[i],errorsVar[i],runtimesL1[i]) = L1Error(s,N[i]);
    println("Error Lasso for N = ",N[i], " is " , errorsL1[i],", runtime is ",runtimesL1[i]);
end

errorsSG = zeros(size(N));
errorsSGL1 = zeros(size(N));
errorsSGExp = zeros(size(N));
errorsSGVar = zeros(size(N));
runtimes = zeros(size(N));

s.filterType = "noFilter";
for i = 1:length(errors)
    (errorsSG[i],errorsSGL1[i],errorsSGExp[i],errorsSGVar[i],runtimes[i]) = L1Error(s,N[i]);
    println("Error SG for N = ",N[i], " is " , errorsSG[i],", runtime is ",runtimes[i]);
end

fig, ax = subplots()
ax[:plot](N,errors, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,errorsSG, "b-+", linewidth=2, label=L"SG", alpha=0.6)
ax[:plot](N,(errorsSG[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("L2 error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsL1, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,errorsSGL1, "b-+", linewidth=2, label=L"SG", alpha=0.6)
ax[:plot](N,(errorsSGL1[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("L1 error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsExp, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,errorsSGExp, "b-+", linewidth=2, label=L"SG", alpha=0.6)
ax[:plot](N,(errorsSGExp[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("Exp error", fontsize=20);

fig, ax = subplots()
ax[:plot](N,errorsVar, "k--o", linewidth=2, label="Lasso Closure", alpha=0.6)
ax[:plot](N,errorsSGVar, "b-+", linewidth=2, label=L"SG", alpha=0.6)
ax[:plot](N,(errorsSGVar[1]*sqrt(NMin))*1.0./sqrt.(N), "r:", linewidth=2, label=L"slope = 1/2", alpha=0.6)
#ax[:plot](N,(errors[1]*NMin)*1.0./N, "g:", linewidth=2, label=L"slope = 1", alpha=0.6)
ax[:legend](loc="lower left")
ax[:set_yscale]("log")
ax[:set_xscale]("log")
ax[:set_xlabel]("N", fontsize=20);
ax[:set_ylabel]("Var error", fontsize=20);

writedlm("results/NValsNx$(s.Nx).txt", N)

writedlm("results/errorsVarNx$(s.Nx).txt", errorsVar)
writedlm("results/errorsSGVarNx$(s.Nx).txt", errorsSGVar)

writedlm("results/errorsExpNx$(s.Nx).txt", errorsExp)
writedlm("results/errorsSGExpNx$(s.Nx).txt", errorsSGExp)

writedlm("results/errorsL1Nx$(s.Nx).txt", errorsL1)
writedlm("results/errorsSGL1Nx$(s.Nx).txt", errorsSGL1)

writedlm("results/errorsNx$(s.Nx).txt", errors)
writedlm("results/errorsSGNx$(s.Nx).txt", errorsSG)

writedlm("results/runtimesSGNx$(s.Nx).txt", runtimes)
writedlm("results/runtimesL1Nx$(s.Nx).txt", runtimesL1)

println("test finished")