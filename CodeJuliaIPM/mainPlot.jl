using NPZ
include("settings.jl")
include("euler.jl")
include("quadrature.jl")
include("Basis.jl")
include("plotting.jl")

close("all");

s = Settings();
q = Quadrature(s.Nq,"Gauss");
basis = Basis(q,s);
Nx = 1000;

x = linspace(s.a,s.b,Nx);

#uSG = npzread("results/solutionFilternoFilterNx700N10tEnd0.0005.jld");
#uL1 = npzread("results/solutionFilterL1FilterClosureNx700N10tEnd0.0005.jld");
uSG = npzread("results/Test2/solutionFilternoFilterNx1000N10tEnd0.14Sigma0.05NewIC.jld");
uL1 = npzread("results/Test2/solutionFilterL1FilterClosureNx1000N10tEnd0.14Sigma0.05NewIC.jld");

tEnd = 0.14;
plotSolution = Plotting(s,basis,q,tEnd);

state = 1;
CompareExpectedValueCons(plotSolution,state,uSG,uL1,x)
ComparePlotInXi(plotSolution,state,uSG,uL1,570)
ComparePlotInXi(plotSolution,state,uSG,uL1,720)
#CompareExpectedValueCons(plotSolution,2,uSG,uL1,x)
#CompareExpectedValueCons(plotSolution,3,uSG,uL1,x)

#uL1 = npzread("results/solutionFilterL1FilterClosureNx1000N15tEnd0.14NewIC.jld");
#PlotExpectedValue(plotSolution,1,uL1)
#PlotExpectedValue(plotSolution,2,uL1)
#PlotExpectedValue(plotSolution,3,uL1)