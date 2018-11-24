using NPZ
include("settings.jl")
include("quadrature.jl")
include("Basis.jl")
include("closure.jl")
include("plotting.jl")


close("all");

s = Settings();
q = Quadrature(s.Nq,"Gauss");
basis = Basis(q,s);
closure = Closure(s,basis,q);
Nx = s.Nx;

x = linspace(s.a,s.b,Nx);

v = npzread("results/Burgers2DNx6000N8tEnd0.01115Sigma0.2.jld");

tEnd = 0.01115;
plotSolution = Plotting(s,basis,q,closure,tEnd);

v0 = 0.5*(s.uL+s.uM+s.sigma0*0.0+s.sigma1*0.0);
xPos = s.x0+v0*tEnd;
xPos = 0.4;
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*xPos)+4)); #0.825 0.797 6385
PlotInXiV(plotSolution,v,convert(Int, round(s.Nx*0.6385)+4));
#PlotExpectedValue(plotSolution,u)