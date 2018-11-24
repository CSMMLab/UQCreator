__precompile__
struct Closure
    uPlus::Float64;
    uMinus::Float64;
    N::Int64;
    Nq::Int64;
    dx::Float64;

    epsNewton::Float64;
    NMaxNewton::Int64;
    basis::Basis;
    q::Quadrature;
    HQuad::Array{Float64,3};
    
    function Closure(settings::Settings,basis::Basis,q::Quadrature)
        Nq = settings.Nq;
        N = settings.N;
        uPlus = settings.uPlus;
        uMinus = settings.uMinus;
        # precompute Hessian terms
        HQuad = zeros(Nq*Nq,N*N,N*N)
        for k = 1:Nq
            for l = 1:Nq
                for i = 1:N
                    for j = 1:N
                        for i1 = 1:N
                            for j1 = 1:N
                                HQuad[(l-1)*Nq+k,(i1-1)*N+i,(j1-1)*N+j] = 0.25*q.w[k]*q.w[l]*basis.PhiQuad[k,i]*basis.PhiQuad[k,j]*basis.PhiQuad[l,i1]*basis.PhiQuad[l,j1];
                            end
                        end
                    end
                end
            end
        end
        new(uPlus,uMinus,N,Nq,settings.dx,settings.epsNewton,settings.NMaxNewton,basis,q,HQuad);
    end
end

function Ubb(obj::Closure,V::Array{Float64,1})
    N = length(V);
    y = zeros(size(V));
    for i = 1:N
        if V[i] > 0
            y[i] = obj.uPlus./(exp(-V[i])+1.0)+obj.uMinus*exp(-V[i])/(1.0+exp(-V[i]));
        else
            y[i] = obj.uMinus./(exp(V[i])+1.0)+obj.uPlus*exp(V[i])/(1.0+exp(V[i]));
        end
    end
    return y;
end

function Ubb(obj::Closure,V::Float64)
    N = length(V);
    y = zeros(size(V));
    if V > 0
        return obj.uPlus./(exp(-V)+1.0)+obj.uMinus*exp(-V)/(1.0+exp(-V));
    else
        return obj.uMinus./(exp(V)+1.0)+obj.uPlus*exp(V)/(1.0+exp(V));
    end
end

function DUbb(obj::Closure, V::Array{Float64,1})
    y = zeros(size(V));
    for k = 1:length(V)
        if V[k] < 0
            y[k] = exp(V[k]).*( obj.uPlus - obj.uMinus )./(1.0+2.0*exp(V[k])+exp(2.0*V[k]));
        else
            y[k] = exp(-V[k]).*( obj.uPlus - obj.uMinus )./(exp(-2.0*V[k])+2.0*exp(-V[k])+1.0);
        end
    end
    return y;
end

function newton(obj::Closure,u, x, nmax, eps)
    xNew = zeros(size(x));
    for n = 1:nmax
        counter::Int64 = 0;
        stepSize::Float64 = 1.0;
        dT = dL(obj,u,x);
        if norm(dT,2) < eps
            return x;
        end
        d = ddL(obj,x)\dT;
        xNew = x - stepSize*d;

        dTNew = dL(obj,u,xNew);
        #println("search direction is ",d)
        #println("residual is ",norm(dT,2));
        while norm(dT,2) < norm(dTNew,2)
            #println("refining...")
            counter = counter+1;
            stepSize = 0.5*stepSize;
            xNew = x - stepSize*d;
            dTNew = dL(obj,u,xNew);
            
            if norm(dTNew,2) < eps
                return xNew;
            end
            if counter > 100 
                error("[ERROR] Newton: Too many stepsize refinements - aborting");
                return xNew;
            end        
        end
        
        if norm(dTNew,2) < eps
            return xNew;
        end
        x = xNew;
    end
    error("[ERROR] Newton: Cannot find optimum!");
    return xNew;
end

function dL(obj::Closure,u::Array{Float64,1},v::Array{Float64,1})
    return ComputeMoments(obj,v) - u;
end

function ComputeMoments(obj::Closure,v::Array{Float64,1})
    return ComputeMoments(obj.basis,0.25*Ubb(obj,EvalAtQuad(obj.basis,v)));
end

function ddL(obj::Closure,v::Array{Float64,1})
    H = zeros(obj.N*obj.N,obj.N*obj.N)
    duVals = DUbb(obj,EvalAtQuad(obj.basis,v));
    for k =1:(obj.Nq*obj.Nq)
        H = H+obj.HQuad[k,:,:]*duVals[k];
    end
    return H;
end

function Solve(obj::Closure, u::Array{Float64,1},v::Array{Float64,1})
    return newton( obj, u, v, obj.NMaxNewton, obj.epsNewton );
end