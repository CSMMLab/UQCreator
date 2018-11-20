include("settings.jl")
include("quadrature.jl")
include("Basis.jl")
include("closure.jl")
include("plotting.jl")

s = Settings();

############################
q = Quadrature(s.Nq,"Gauss"); # ClenshawCurtis Gauss
b = Basis(q,s)
closure = Closure(s,b,q);

function Filter(u::Array{Float64,1})
    uF = zeros(size(u));
    if s.filterType == "L2Filter"
        for i = 1:s.N
            uF[i] = u[i]/(1.0+s.lambdaFilter*i^2*(i-1)^2);
        end
    elseif s.filterType == "L1Filter"
        for i = 1:s.N
            scL1 = 1.0-4.0*s.lambdaFilter*i*(i-1)*b.PhiL1[i]/abs(u[i]); # 4.0 term to scale filter strength, can be left out
            if scL1 < 0
                scL1 = 0;
            end
            uF[i] = scL1*u[i];
        end
    elseif s.filterType == "L1FilterClosure"
        N = s.N;
        lambdaFilter = abs(u[N])/(N*(N-1)*b.PhiL1[N]);
        for i = 2:s.N
            scL1 = 1.0-lambdaFilter*i*(i-1)*b.PhiL1[i]/abs(u[i]);
            if scL1 < 0 || abs(u[i]) < 1.0e-7
                scL1 = 0;
            end
            uF[i] = scL1*u[i];
        end
    end
    return uF;
end

function IC(x,a,b)
    #return sin(x*pi/s.dx)+7.5;
    y = zeros(size(x));
    deltaU = 0.00001;
    uL = 12.0-deltaU;
    uR = 1.0+deltaU;
    dx = b-a;
    x0 = a+0.3*dx;
    x1 = x0;
    #x1 = a+0.6*dx;
    for k = 1:length(x) 
        if x[k] < x0
            y[k] = uL;
        elseif x[k] < x1
            y[k] = uL + (uR - uL)*(x[k]-x0)/(x1-x0);
        else
            y[k] = uR;
        end
    end
    return y;
end

j = 4;

u = zeros(s.N);
for i = 1:s.N
    u[i] = Integral(q, xVal->IC(xVal,-1.0,1.0).*Phi(i-1,xVal)*0.5,-1.0,1.0);
end


uF = Filter(u);

lambda = zeros(s.N);
y = Solve(closure,uF,lambda)
xFine = linspace(-1.0,1.0,100);
xCoarse = q.xi;

fig, ax = subplots()
ax[:plot](xFine,IC(xFine,-1.0,1.0), "k-", linewidth=2, label=L"$u_{exact}$", alpha=0.6)
ax[:plot](xCoarse, Ubb(closure,EvalAtQuad(b,y)), "r--", linewidth=2, label=L"$u_{IPM}$", alpha=0.6)
ax[:plot](xCoarse,EvalAtQuad(b,u), "b:", linewidth=2, label=L"$u_{SG}$", alpha=0.6)
ax[:plot](xCoarse,EvalAtQuad(b,uF), "g-", linewidth=2, label=L"$u_{fSG}$", alpha=0.6)
ax[:legend](loc="upper right")

uTest = zeros(s.N);
for i = 1:s.N
    uTest[i] = IntegralVec(q, 0.5*Ubb(closure,EvalAtQuad(b,y)).*b.PhiQuad[:,i],-1.0,1.0);
end

println("Moment Vector Test ",uTest);
println("Moment Vector ",u);

println("Sol at quad! -> ", Ubb(closure,EvalAtQuad(b,y)));
println("DONE! -> ", y);
############################