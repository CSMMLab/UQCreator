__precompile__
import GSL

struct Basis
    # number of moments
    N::Int64;
    NXi::Int64;
    Nq::Int64;
    
    # precomputed Legendre Polynomials at quadrature points
    PhiQuad::Array{Float64,2};
    PhiQuadW::Array{Float64,2};
    PhiL1::Array{Float64,1}

    function Basis(quadrature::Quadrature,settings::Settings)
        N = settings.N;
        Nq = settings.Nq;
        NXi = settings.NXi;

        # precompute Legendre basis functions at quad points
        PhiQuad = zeros(Nq*Nq,N*N);
        PhiQuadW = zeros(N*N,Nq*Nq);
        PhiTQuad = zeros(N*N,Nq*Nq);
        PhiL1 = zeros(N*N);
        for i = 1:N
            for j = 1:N
                #PhiL1[i] = Integral(quadrature,xi->abs(Phi(i-1,xi))*0.5,-1.0,1.0 ); # TODO: Add for L1 filter
                for k = 1:Nq
                    for q = 1:Nq
                        PhiQuad[(q-1)*Nq+k,(j-1)*N+i] = Phi.(i-1,quadrature.xi[k])*Phi.(j-1,quadrature.xi[q]);
                        PhiTQuad[(j-1)*N+i,(q-1)*Nq+k] = Phi.(i-1,quadrature.xi[k])*Phi.(j-1,quadrature.xi[q]);
                        PhiQuadW[(j-1)*N+i,(q-1)*Nq+k] = PhiTQuad[(j-1)*N+i,(q-1)*Nq+k]*quadrature.w[k]*quadrature.w[q];
                    end
                end
            end
        end
        for k = 1:(Nq^2)
            println(PhiQuad[k,:])
        end
        assert(false)

        new(N,Nq,NXi,PhiQuad,PhiQuadW,PhiL1);
    end

end

# Legendre Polynomials on [-1,1]
function Phi(n::Int64,xi)
    return sqrt(2.0*n+1.0)*GSL.sf_legendre_Pl.(n,xi);;
end

# MultiD Legendre Polynomials on [-1,1]
function Phi(obj::Basis,n::Array{Int64,1},xi::Array{Float64,1})
    tmp = 1.0;
    for i = 1:obj.NXi
        tmp = tmp * sqrt(2.0*n[i]+1.0)*GSL.sf_legendre_Pl.(n[i],xi[i]);
    end
    return tmp;
end

# evaluate polynomial with moments u at spatial position x
function Eval(obj::Basis,u::Array{Float64,1},xi,eta)
    assert(length(xi) == length(eta))
    y=zeros(length(xi))
    for i = 1:obj.N
        for j = 1:obj.N
            y = y+u[(j-1)*obj.N+i]*Phi.(i-1,xi).*Phi.(j-1,eta);
        end
    end
    return y;
end

function EvalAtQuad(obj::Basis,u::Array{Float64,1})
    return obj.PhiQuad*u;
end

# returns N moments of all states evaluated at the quadrature points and stored in uQ
# uQ in [Nq,states], returns [N,states]
function ComputeMoments(obj::Basis,uQ::Array{Float64,1})
    return obj.PhiQuadW*uQ;
end