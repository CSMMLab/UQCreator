clc
clearvars
close all

problem = 2;                                    % 1 = Euler ideal gas; 2 = Burgers; 3 = Advection; 0 is no longer in the code ;)
bounds = [0 10];
nCells = 100;
grid = linspace(bounds(1), bounds(2),nCells)';
ic = zeros(nCells,1);                           % initial condition
ic(grid<5) = 1;                                 % small shock
dt = 0.1;
tEnd = 2;
nTimeSteps = tEnd/dt;
maxCFL = 0.5;
BC = 0;                                         % 0 = open; 1 = periodic
limiter = 1;                                    % 0 = none; 1 = vanLeer; 2 = MinMod
 
tic
res = FVSolver1D(problem, ic, dt, nTimeSteps, maxCFL, BC, limiter, bounds);
time = toc

plot(res(:,1), res(:,2));
