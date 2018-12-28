function IC1Exact(t,x,xi,sigma,uL,uR)
    x0Save = 0.5;
    x1Save = 1.5;

    if t >= (x1Save-x0Save)/(uL-uR);
        tS = (x1Save-x0Save)/(uL-uR);
        x0BeforeShock = x0Save+sigma*xi + tS*uL;
        x1BeforeShock = x1Save+sigma*xi + tS*uR;
        x0 = x0BeforeShock + (t-tS)*(uL+uR)*0.5;
        x1 = x0 - 1.0;
    else
        x0 = x0Save+sigma*xi + t*uL;
        x1 = x1Save+sigma*xi + t*uR;
    end
    if x < x0
        return uL;
    elseif x < x1
        return uL + (uR - uL)*(x-x0)/(x1-x0);
    else
        return uR;
    end
    
end

using PyPlot

Nx = 1000;
a = 0;
b = 3;
x = linspace(a,b,Nx);
sigmaRand = 1.0;
sigma = 0.0;
uL = 12.0;
uR = 3.0;
mu = 0;
NSamples = 1000;
tmp = zeros(Nx);
t = 0.12;

xi = 1
IC1Exact(t,x[2],xi,sigma,uL,uR)

for k = 1:NSamples
    xi = rand(); # Gauss Distr
    #xi = 2*rand()-1;
    for j = 1:Nx
        tmp[j] = tmp[j] + IC1Exact(t,x[j],xi,sigma,uL,uR)
    end
end

tmp = tmp/NSamples;

# plot
fig, ax = subplots(figsize=(15, 8), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax[:plot](x,tmp, "r--", linewidth=2, label=L"$E[u_{N}]$", alpha=0.6)
ylimMinus = -0.5;
ylimPlus = 16.0
#ax[:set_ylim]([ylimMinus,ylimPlus])
ax[:set_xlim]([a,b])
ax[:set_xlabel]("x", fontsize=20);
