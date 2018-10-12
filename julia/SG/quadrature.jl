__precompile__
import FastTransforms

struct Quadrature
    Nq::Int;
    w::Array{Float64,1};
    xi::Array{Float64,1};

   function Quadrature(NqVal,quadType)
        if quadType == "Gauss"
            xi,w = QuadNodesWeights(NqVal,-1.0,1.0)
        elseif quadType == "ClenshawCurtis"
            xi,w = QuadNodesWeightsClCu(NqVal,-1.0,1.0)
        end
        xi = xi[end:-1:1];
        w = w[end:-1:1];

        new(NqVal,w,xi)
   end

end

function QuadNodesWeightsClCu(Nq::Int64,a::Float64,b::Float64) # transformation to a,b not tested!
    xiQuad,wQuad = FastTransforms.clenshawcurtis(Nq,0.,0.)
    xiQuad=(a*(1-xiQuad)+b*(1+xiQuad))/2.0;
    wQuad = wQuad*(2.0/(b-a));
    return xiQuad,wQuad;
end

# compute Gauss nodes and weights
function QuadNodesWeights(N::Int64,a::Float64,b::Float64)

    N=N-1;
    N1=N+1; N2=N+2;

    xu=linspace(-1,1,N1);

    # Initial guess
    y=cos.((2*(0:N)+1)*pi/(2*N+2))+(0.27/N1)*sin.(pi*xu*N/N2);

    # Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);

    # Derivative of LGVM
    Lp=zeros(N1,N2);

    # Compute the zeros of the N+1 Legendre Polynomial
    # using the recursion relation and the Newton-Raphson method

    y0=2;

    eps = 2.2204e-16;

    # Iterate until new points are uniformly within epsilon of old points
    while maximum(abs.(y-y0)) > eps
        L[:,1]=1;
        Lp[:,1]=0;

        L[:,2]=y;
        #Lp[:,2]=1;

        for k=2:N1
            L[:,k+1]=( (2*k-1)*y.*L[:,k]-(k-1)*L[:,k-1] )/k;
        end

        Lp=(N2)*( L[:,N1]-y.*L[:,N2] )./(1-y.^2);

        y0=y;
        y=y0-(L[:,N2]./Lp);

    end

    # Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;

    # Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

    return x,w;
end

function Integral(q::Quadrature, f)
    dot( q.w, f.(q.xi) );
end

function Integral(q::Quadrature, f, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    xi = (a*(1.0-q.xi)+b*(1+q.xi))/2.0;
    #return q.w'f.(xi);
    return sum(w'f.(xi))
end

function IntegralVec(q::Quadrature, fVec::Array{Float64,1}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return sum(w'fVec)
end

function IntegralVec(q::Quadrature, fVec::Array{Float64,1})
    return sum(q.w'fVec)
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2}, a::Float64, b::Float64)
    w = q.w*(b-a)*0.5;
    return fVec*w;
end

function IntegralMat(q::Quadrature, fVec::Array{Float64,2})
    return fVec*q.w;
end