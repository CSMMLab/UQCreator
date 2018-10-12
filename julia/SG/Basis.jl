__precompile__
import GSL

struct Basis
    # number of moments
    N::Int64;
    Nq::Int64;

    # precomputed Legendre Polynomials at quadrature points
    PhiQuad::Array{Float64,2};
    PhiTQuad::Array{Float64,2};
    PhiQuadW::Array{Float64,2};
    PhiL1::Array{Float64,1}

    function Basis(quadrature::Quadrature,settings::Settings)
        N = settings.N;
        Nq = settings.Nq;

        # precompute Legendre basis functions at quad points
        PhiQuad = zeros(settings.Nq,N);
        PhiQuadW = zeros(N,settings.Nq);
        PhiTQuad = zeros(N,settings.Nq);
        PhiL1 = zeros(N);
        for i = 1:N
            PhiL1[i] = Integral(quadrature,xi->abs(Phi(i-1,xi))*0.5,-1.0,1.0 );
            for k = 1:settings.Nq
                PhiQuad[k,i] = Phi.(i-1,quadrature.xi[k]); 
                PhiTQuad[i,k] = Phi.(i-1,quadrature.xi[k]);
                PhiQuadW[i,k] = PhiTQuad[i,k]*quadrature.w[k];
            end
        end

        new(N,Nq,PhiQuad,PhiTQuad,PhiQuadW,PhiL1);
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

# evalueates states at quadrature points. returns solution[quadpoints,states]
function EvalAtQuad(obj::Basis,u::Array{Float64,2})
    return obj.PhiQuad*u;
end

# returns N moments of a function evaluated at the quadrature points and stored in uQ
function ComputeMoments(obj::Basis,uQ::Array{Float64,1})
    return obj.PhiQuadW*uQ;
end

# returns N moments of all states evaluated at the quadrature points and stored in uQ
# uQ in [Nq,states], returns [N,states]
function ComputeMoments(obj::Basis,uQ::Array{Float64,2})
    return obj.PhiQuadW*uQ;
end