include("settings.jl")
include("quadrature.jl")
include("Basis.jl")
include("closure.jl")
include("plotting.jl")

s = Settings();

############################
q = Quadrature(s.Nq,"Gauss");
b = Basis(q,s);
closure = Closure(s,b,q);

function IC(xi::Array{Float64,1})
    return sin(xi[1]*pi) + exp(0.1*xi[2]);
end

function mymod(k,N)
    return mod(k-1,N)+1;
end

function Indices(k,N,dim)
    ind::Array{Int,1} = zeros(dim)
    for l = 0:(dim-1)
        ind[l+1] = mod( Int((k-mymod(k,s.Nq^l))/(s.Nq^l)),s.Nq)+1;
    end
    return ind;
end

# setup quadpoints
xiQ = zeros(2,s.Nq*s.Nq);
xiQ2 = zeros(2,s.Nq*s.Nq);
xi = q.xi;

for l = 1:s.Nq
    for k = 1:s.Nq
        xiQ[:,(l-1)*s.Nq+k] = [xi[l],xi[k]];
    end
end

dim = 4;
for k = 1:(s.Nq^dim)
    i = Indices(k,s.Nq,dim);
    modVal1 = mymod(k,s.Nq);
    modVal2 = mod( Int((k-mymod(k,s.Nq^1))/(s.Nq^1)),s.Nq)+1;
    modVal3 = mod( Int((k-mymod(k,s.Nq^2))/(s.Nq^2)),s.Nq)+1;
    modVal4 = mod( Int((k-mymod(k,s.Nq^3))/(s.Nq^3)),s.Nq)+1;
    println("mod = ",modVal4," ",i[4],", ",modVal3," ",i[3],", ",modVal2," ",i[2],", ",modVal1," ",i[1])
    #xiQ2[:,k] = [xi[Int((k-modVal)/s.Nq) + 1],xi[modVal]];
end

println("distance is ",maximum(xiQ2-xiQ))

A = zeros(s.N*s.N,s.Nq*s.Nq)
for k = 1:s.Nq
    for l = 1:s.Nq
        for j = 1:s.N
            for i = 1:s.N
                A[(j-1)*s.N+i,(l-1)*s.Nq+k] = b.PhiQuad[l,i]*b.PhiQuad[k,j]*q.w[k]*q.w[l];
            end
        end
    end
end

IC(xiQ[:,2])