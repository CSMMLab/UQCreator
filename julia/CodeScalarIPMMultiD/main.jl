using Revise
using NPZ
include("settings.jl")
include("Solver.jl")
include("plotting.jl")

close("all")

#pointer_from_objref(x);

s = Settings();

s.closureType = "BoundedBarrier"

############################
solver = IPMSolver(s);
#solver = SolverOnMoments(s);

@time tEnd, u, v = Solve(solver);

plotSolution = Plotting(s,solver.basis,solver.q,solver.closure,tEnd);

#PlotInX(plotSolution,u,[-1.0, 0.0, 1.0]);
v0 = 0.5*(s.uL+s.uM+s.sigma0*0.0+s.sigma1*0.0);
xPos = s.x0+v0*tEnd;
xPos = 0.4;
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*xPos)+4)); #0.825 0.797 6385
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*0.6385)+4));
#PlotInXi(plotSolution,u,convert(Int, round(s.Nx*1.0/3.0)));
PlotExpectedValue(plotSolution,u)

npzwrite("results/Burgers2DClosure$(s.closureType)Nx$(s.Nx)N$(s.N)Nq$(s.Nq)tEnd$(s.tEnd)Sigma$(s.sigma0).jld", v)

s.closureType = "LogBarrier"
s.uMinus = 0.0;
s.uPlus = 12.5;

############################
solver = IPMSolver(s);
#solver = SolverOnMoments(s);

@time tEnd, u, v = Solve(solver);

plotSolution = Plotting(s,solver.basis,solver.q,solver.closure,tEnd);

#PlotInX(plotSolution,u,[-1.0, 0.0, 1.0]);
v0 = 0.5*(s.uL+s.uM+s.sigma0*0.0+s.sigma1*0.0);
xPos = s.x0+v0*tEnd;
xPos = 0.4;
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*xPos)+4)); #0.825 0.797 6385
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*0.6385)+4));
#PlotInXi(plotSolution,u,convert(Int, round(s.Nx*1.0/3.0)));
PlotExpectedValue(plotSolution,u)

npzwrite("results/Burgers2DClosure$(s.closureType)Nx$(s.Nx)N$(s.N)Nq$(s.Nq)tEnd$(s.tEnd)Sigma$(s.sigma0)DeltaUPoette.jld", v)

s.uMinus = 1.0-0.001;
s.uPlus = 12.2+0.001;

############################
solver = IPMSolver(s);
#solver = SolverOnMoments(s);

@time tEnd, u, v = Solve(solver);

plotSolution = Plotting(s,solver.basis,solver.q,solver.closure,tEnd);

#PlotInX(plotSolution,u,[-1.0, 0.0, 1.0]);
v0 = 0.5*(s.uL+s.uM+s.sigma0*0.0+s.sigma1*0.0);
xPos = s.x0+v0*tEnd;
xPos = 0.4;
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*xPos)+4)); #0.825 0.797 6385
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*0.6385)+4));
#PlotInXi(plotSolution,u,convert(Int, round(s.Nx*1.0/3.0)));
PlotExpectedValue(plotSolution,u)

npzwrite("results/Burgers2DClosure$(s.closureType)Nx$(s.Nx)N$(s.N)Nq$(s.Nq)tEnd$(s.tEnd)Sigma$(s.sigma0)DeltaUSmall.jld", v)

println("test finished")
