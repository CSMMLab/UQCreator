% run mex -setup if you haven't already...

SRC_DIR = fullfile('..','matlabFVSolver1D');
INCLUDES = SRC_DIR;
SRC = {...    
    'Cell.cpp'...
    'Cell1D.cpp'...
    'Advection.cpp'...
    'BurgersEquation.cpp'...
    'EulerEquationsIdGas.cpp'...
    'Heun.cpp'...
    'HLLSolver1D.cpp'...
    'Limiter.cpp'...
    'Mesh.cpp'...
    'Mesh1D.cpp'...
    'MinMod1D.cpp'...
    'PhysicalProblem.cpp'...
    'Rhs.cpp'...
    'Rhs1D.cpp'...
    'RiemannSolver.cpp'...
    'Settings.cpp'...
    'TimeIntegrator.cpp'...
    'VanLeer1D.cpp'...
    'Vector.cpp'};
for i = 1:length(SRC)
    SRC{i} = fullfile(SRC_DIR,SRC{i});
end

mex('FVSolver1D.cpp', '-I/usr/include', '-UDEBUG', ['-I',SRC_DIR], SRC{:} );