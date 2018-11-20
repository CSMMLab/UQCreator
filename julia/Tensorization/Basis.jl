__precompile__
import GSL

struct Basis
    # number of moments
    N::Int64;
    Nq::Int64;

    # precomputed Legendre Polynomials at quadrature points
    PhiQuad::Array{Float64,2};
    PhiL1::Array{Float64,1}

    function Basis(quadrature::Quadrature,settings::Settings)
        N = settings.N;
        Nq = quadrature.Nq;

        # precompute Legendre basis functions at quad points
        PhiQuad = zeros(Nq,N);
        PhiL1 = zeros(N);
        for i = 1:N
            PhiL1[i] = Integral(quadrature,xi->abs(Phi(i-1,xi))*0.5,-1.0,1.0 );
            for k = 1:Nq
                PhiQuad[k,i] = Phi.(i-1,quadrature.xi[k]);
            end
        end

        new(N,Nq,PhiQuad,PhiL1);
    end
end

# Legendre Polynomials on [-1,1]
function Phi(n::Int64,xi::Float64) ## can I use this with input vectors xVal?
    return sqrt(2.0*n+1.0)*GSL.sf_legendre_Pl.(n,xi);
end

# evaluate polynomial with moments u at spatial position x
function Eval(obj::Basis,u::Array{Float64,1},xi)
    y=zeros(length(xi))
    for i = 1:obj.N
        y = y+u[i]*Phi.(i-1,xi);
    end
    return y;
end

function EvalAtQuad(obj::Basis,u::Array{Float64,1})
    return obj.PhiQuad*u;
end

# functions for evaluating entropic variables

# evaluate dual State at quadrature points for given dual variables v
function EvalVAtQuad(obj::Basis,v::Array{Float64,1})
    return obj.PhiTildeQuad*v;
end