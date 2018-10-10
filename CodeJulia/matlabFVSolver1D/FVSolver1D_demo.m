clc
clearvars
close all

rhoL = 1.0; rhoR = 0.5; %0.1;
pL = 1.0; pR = 0.65;
uL = 0.0; uR = 0.0;
Nx = 500;
x0 = 1.0;

problem = 1;                                    % 1 = Euler ideal gas; 2 = Burgers; 3 = Advection; 0 is no longer in the code ;)
bounds = [0 3];
grid = linspace(bounds(1), bounds(2),Nx)';
ic = zeros(Nx,3);                           % initial condition
ic(grid<x0,1) = rhoL;                                 % small shock
ic(grid>=x0,1) = rhoR;                                 % small shock
ic(grid<x0,2) = pL; 
ic(grid>=x0,2) = pR;
ic(grid<x0,3) = uL; 
ic(grid>=x0,3) = uR;
dt = 1.2/Nx;
tEnd = 0.5;
nTimeSteps = tEnd/dt;
maxCFL = 0.5;
BC = 0;                                         % 0 = open; 1 = periodic
limiter = 0;                                    % 0 = none; 1 = vanLeer; 2 = MinMod
  
res = FVSolver1D(problem, ic, dt, nTimeSteps, maxCFL, BC, limiter, bounds);

size(res)
plot(res(:,1), res(:,2),res(:,1), res(:,3),res(:,1), res(:,4));
