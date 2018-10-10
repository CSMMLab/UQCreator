__precompile__
struct Closure
    N::Int64;
    Nq::Int64;
    dx::Float64;
    states::Int64;
    gamma::Float64;

    epsNewton::Float64;
    NMaxNewton::Int64;
    basis::Basis;
    q::Quadrature;
    HQuad::Array{Float64,3};
    
    function Closure(settings::Settings,basis::Basis,q::Quadrature)

        # precompute Hessian terms
        HQuad = zeros(q.Nq,basis.N,basis.N)
        for k = 1:q.Nq
            for i = 1:basis.N
                for j = 1:basis.N
                    HQuad[k,i,j] = 0.5*q.w[k]*basis.PhiQuad[k,i]*basis.PhiQuad[k,j];
                end
            end
        end
        new(basis.N,q.Nq,settings.dx,settings.states,settings.gamma,settings.epsNewton,settings.NMaxNewton,basis,q,HQuad);
    end
end

function UKin(obj::Closure,V::Array{Float64,2})
    y = zeros(size(V));
    gamma = obj.gamma;
    for k = 1:obj.Nq
        v1 = V[k,1]; v2 = V[k,2]; v3 = V[k,3]; v4 = V[k,4];
        expTerm = (-exp(((v2^2 + v3^2 - 2*v1*v4 - 2*v4*gamma)/(2*v4)))*v4)^(1/(1 + gamma));
        y[k,1] = expTerm;
        y[k,2] = -((v2*expTerm)/v4);
        y[k,3] = -((v3*expTerm)/v4);
        y[k,4] = -((expTerm*(-v2^2 - v3^2 + 2*v4))/(2*v4^2));
    end
    return y;
end

function UKin(obj::Closure,V::Array{Float64,1})
    y = zeros(3);
    gamma = obj.gamma;
    v1 = V[1]; v2 = V[2]; v3 = V[3]; v4 = V[4];
    expTerm = (-exp(((v2^2 + v3^2 - 2*v1*v4 - 2*v4*gamma)/(2*v4)))*v4)^(1/(1 + gamma));
    y[1] = expTerm;
    y[2] = -((v2*expTerm)/v4);
    y[3] = -((v3*expTerm)/v4);
    y[4] = -((expTerm*(-v2^2 - v3^2 + 2*v4))/(2*v4^2));

    return y;
end

function DUKin(obj::Closure, V::Array{Float64,2})
    y = zeros(obj.Nq,obj.states,obj.states);
    gamma = obj.gamma;
    for k = 1:obj.Nq
        v1 = V[k,1]; v2 = V[k,2]; v3 = V[k,3]; v4 = V[k,4];
        expTerm = (-exp(((v2^2 + v3^2 - 2*v4*(v1 + gamma))/(2*v4)))*v4)^(1/(1 + gamma))*1/(1 + gamma);
        vTerm = v2^2 + v3^2 + 2*v4*gamma;
        y[k,1,1] = -(expTerm);
        y[k,1,2] = (v2*expTerm)/(v4);
        y[k,1,3] = (v3*expTerm)/(v4);
        y[k,1,4] = -(((v2^2 + v3^2 - 2*v4)*expTerm)/(2*v4^2));
        y[k,2,1] = (v2*expTerm)/(v4);
        y[k,2,2] =  -((expTerm*(v2^2 + v4 + v4*gamma))/(v4^2));
        y[k,2,3] =  -((v2*v3*expTerm)/(v4^2));
        y[k,2,4] =  (v2*expTerm*vTerm)/(2*v4^3);
        y[k,3,1] =  (v3*expTerm)/(v4);
        y[k,3,2] =  -((v2*v3*expTerm)/(v4^2)); 
        y[k,3,3] =  -((expTerm*(v3^2 + v4 + v4*gamma))/(v4^2)); 
        y[k,3,4] =  (v3*expTerm *vTerm)/(2*v4^3 );
        y[k,4,1] =  -(((v2^2 + v3^2 - 2*v4)*expTerm)/(2*v4^2));
        y[k,4,2] =  (v2*expTerm *vTerm)/(2*v4^3 );
        y[k,4,3] =  (v3*expTerm *vTerm)/(2*v4^3); 
        y[k,4,4] =  -(expTerm*(v2^4 + v3^4 + 4*v3^2*v4*gamma - 4*v4^2*gamma + 2*v2^2*(v3^2 + 2*v4*gamma)))/(4*v4^4);
    end
    return y;
end

function MakeVector(obj::Closure,matrix::Array{Float64,2})
    y = zeros(obj.states*obj.N)
    for s = 1:obj.states
        for i = 1:obj.N
            y[(s-1)*obj.N+i] = matrix[i,s];
        end
    end
    return y;
end

function MakeMatrix(obj::Closure,vector::Array{Float64,1})
    y = zeros(obj.N,obj.states)
    for s = 1:obj.states
        for i = 1:obj.N
            y[i,s] = vector[(s-1)*obj.N+i];
        end
    end
    return y;
end

function newton(obj::Closure,u, x, nmax, eps)
    xNew = zeros(size(x));
    for n = 1:nmax
        counter::Int64 = 0;
        stepSize::Float64 = 0.1;
        dT = dL(obj,u,x);
        if norm(dT,2) < eps
            return x;
        end
        d = ddL(obj,x)\dT;
        dMat = MakeMatrix(obj,d);
        xNew = x - stepSize*dMat;
        #println("moments are ",ComputeMoments(obj,x));
        #println("but should be ",u);
        #println("gradient is ",dT);
        #println("Newton direction is ",dMat);
        dTNew = dL(obj,u,xNew);
        #println("search direction is ",d)
        #println("residual is ",norm(dT,2));
        while norm(dT,2) < norm(dTNew,2)
            #println("refining...")
            counter = counter+1;
            stepSize = 0.5*stepSize;
            xNew = x - stepSize*dMat;
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

function dL(obj::Closure,u::Array{Float64,2},v::Array{Float64,2})
    y = zeros(obj.N,obj.states);
    UKinVals = UKin(obj,EvalAtQuad(obj.basis,v));
    for s = 1:obj.states
        y[:,s] = IntegralMat(obj.q,0.5*(obj.basis.PhiQuad').*UKinVals[:,s]') - u[:,s];
    end
    return MakeVector(obj,y);
end

function ComputeMoments(obj::Closure,v::Array{Float64,2})
    y = zeros(obj.N,obj.states);
    UKinVals = UKin(obj,EvalAtQuad(obj.basis,v));
    for s = 1:obj.states
        y[:,s] = IntegralMat(obj.q,0.5*(obj.basis.PhiQuad').*UKinVals[:,s]');
    end
    return y;
end

function ddL(obj::Closure,v::Array{Float64,2})
    H = zeros(obj.states*obj.N,obj.states*obj.N)
    duVals = DUKin(obj,EvalAtQuad(obj.basis,v));

    for j = 1:obj.N
        for i = 1:obj.N
            for l = 1:obj.states
                for s = 1:obj.states
                    for k = 1:obj.Nq
                        H[(s-1)*obj.N+i,(l-1)*obj.N+j] = H[(s-1)*obj.N+i,(l-1)*obj.N+j] + obj.HQuad[k,i,j]*duVals[k,s,l];
                    end
                end
            end
        end
    end

    return H;
end

function Solve(obj::Closure, u::Array{Float64,2},v::Array{Float64,2})
    #println("u = ",u)
    return newton( obj, u, v, obj.NMaxNewton, obj.epsNewton );
end