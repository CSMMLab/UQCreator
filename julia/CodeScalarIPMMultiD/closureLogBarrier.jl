__precompile__
struct ClosureLogBarrier
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
    
    function ClosureLogBarrier(settings::Settings,basis::Basis,q::Quadrature)
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

function Ubb(obj::ClosureLogBarrier,V::Array{Float64,1})
    N = length(V);
    y = zeros(size(V));
    for i = 1:N
        y[i] = 0.5*(obj.uMinus+obj.uPlus+(obj.uMinus-obj.uPlus)^2*V[i]/(sqrt((obj.uMinus-obj.uPlus)^2*V[i]^2+4)+2));
    end
    return y;
end

function Ubb(obj::ClosureLogBarrier,V::Float64)
    return 0.5*(obj.uMinus+obj.uPlus+(obj.uMinus-obj.uPlus)^2*V/(sqrt((obj.uMinus-obj.uPlus)^2*V^2+4)+2));
end

function dds2(obj::ClosureLogBarrier,u)
    y = 1.0./((u-obj.uMinus).^2)+1.0./((obj.uPlus-u).^2);
end

function DUbb(obj::ClosureLogBarrier, V::Array{Float64,1})
    return 1.0./(dds2(obj,Ubb(obj,V)));
end

function newton(obj::ClosureLogBarrier,u, x, nmax, eps)
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

function dL(obj::ClosureLogBarrier,u::Array{Float64,1},v::Array{Float64,1})
    return ComputeMoments(obj,v) - u;
end

function ComputeMoments(obj::ClosureLogBarrier,v::Array{Float64,1})
    return ComputeMoments(obj.basis,0.25*Ubb(obj,EvalAtQuad(obj.basis,v)));
end

function ddL(obj::ClosureLogBarrier,v::Array{Float64,1})
    H = zeros(obj.N*obj.N,obj.N*obj.N)
    duVals = DUbb(obj,EvalAtQuad(obj.basis,v));
    for k =1:(obj.Nq*obj.Nq)
        H = H+obj.HQuad[k,:,:]*duVals[k];
    end
    return H;
end

function Solve(obj::ClosureLogBarrier, u::Array{Float64,1},v::Array{Float64,1})
    return newton( obj, u, v, obj.NMaxNewton, obj.epsNewton );
end