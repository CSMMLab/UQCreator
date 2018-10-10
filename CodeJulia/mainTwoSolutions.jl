using NPZ
include("settings.jl")
include("mesh.jl")
include("Solver.jl")
include("plotting.jl")


#pointer_from_objref(x);
#close("all");

s = Settings();
m = Mesh(s);

s.filterType = "L1FilterClosure"

############################

solver = Solver(s);

@time tEnd, u = Solve(solver);

npzwrite("results/solutionDuctFilter$(s.filterType)Nx$(s.Nx)N$(s.N)tEnd$(s.tEnd)Sigma$(s.sigma).jld", u)

s.filterType = "noFilter"

solverSG = Solver(s);
@time tEnd, u = Solve(solverSG);

npzwrite("results/solutionDuctFilter$(s.filterType)Nx$(s.Nx)N$(s.N)tEnd$(s.tEnd)Sigma$(s.sigma).jld", u)

println("main finished")
