using NPZ
include("settings.jl")
include("mesh.jl")
include("Solver.jl")
include("plotting.jl")


s = Settings();
m = Mesh(s);

############################

solver = Solver(s);

@time tEnd, u = Solve(solver);

npzwrite("results/solutionDuctFilter$(s.filterType)Nx$(s.Nx)N$(s.N)tEnd$(s.tEnd)Sigma$(s.sigma).jld", u)
