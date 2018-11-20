using Revise
using NPZ
include("settings.jl")
include("IPMSolver.jl")
#include("closure.jl")
include("plotting.jl")

#pointer_from_objref(x);

close()

s = Settings();

############################
solver = IPMSolver(s);

@time tEnd, u, v = Solve(solver);

plotSolution = Plotting(s,solver.basis,solver.q,solver.closure,tEnd);

PlotInX(plotSolution,v,[-1.0, 0.0, 1.0]);
PlotInXi(plotSolution,v,convert(Int, round(s.Nx*0.57)));
#PlotInXi(plotSolution,v,1);
PlotExpectedValue(plotSolution,v)

l1err = L1Error(plotSolution,v,tEnd);
l2err = L2Error(plotSolution,v,tEnd);
errorExp,errorVar = L1ErrorExpVar(plotSolution,v,tEnd);
println("Error L1 IPM for N = ",s.N, " is " , l1err);
println("Error L2 IPM for N = ",s.N, " is " , l2err);
println("Error EXP IPM for N = ",s.N, " is " , errorExp);
println("Error Var IPM for N = ",s.N, " is " , errorVar);

npzwrite("results/BurgersIPMNx$(s.Nx)N$(s.N)tEnd$(s.tEnd)Sigma$(s.sigma).jld", v)

println("test finished")