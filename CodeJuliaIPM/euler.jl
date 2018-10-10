__precompile__

mutable struct Euler
    gamma::Float64; # isotropic exponent
    SpecificR::Float64; # specific gas constant
    function Euler(s::Settings)
        new(s.gamma,s.SpecificR)
    end
end

# determine primitives from conservatives
function primFromCons(obj::Euler,conservatives::Array{Float64,1})
    rho=conservatives[1];
    rhoInv=1/rho;
    u=conservatives[2]*rhoInv;
    v=conservatives[3]*rhoInv;
    kineticEnergy=0.5*rho*(u^2+v^2);
    innerEnergy=(conservatives[4]-kineticEnergy)*rhoInv;
    p=rho*(obj.gamma-1)*innerEnergy;
    a=sqrt(obj.gamma*p*rhoInv);
    T=p*rhoInv/obj.SpecificR;
    return rho,p,u,v,a,T;
end

# determine primitives from conservatives
function SpeedFromCons(obj::Euler,conservatives::Array{Float64,1})
    rhoInv=1/conservatives[1];
    u=conservatives[2]*rhoInv;
    v=conservatives[3]*rhoInv;
    p=(obj.gamma-1.0)*(conservatives[4]-0.5*conservatives[1]*(u^2+v^2));
    a=sqrt(obj.gamma*p*rhoInv);
    return u,v,a;
end

# determine primitives from conservatives
function SpeedPressureFromCons(obj::Euler,conservatives::Array{Float64,1})
    rhoInv=1/conservatives[1];
    u=conservatives[2]*rhoInv;
    v=conservatives[3]*rhoInv;
    p=(obj.gamma-1.0)*(conservatives[4]-0.5*conservatives[1]*(u^2+v^2));
    return p,u,v;
end

# determine conservatives from primitives
function consFromPrim(obj::Euler,primitives::Array{Float64,1})
    rho = primitives[1]; p = primitives[2]; u = primitives[3]; v = primitives[4]
    rhou=rho*u;
    rhov=rho*v;
    kineticEnergy=0.5*rho*(u^2+v^2);
    innerEnergy=(p/(rho*(obj.gamma-1)))*rho;
    rhoe=kineticEnergy+innerEnergy;
    conservatives = [rho;rhou;rhov;rhoe];
    return conservatives;
end

# physical flux
function physicalFlux(obj::Euler,conservatives::Array{Float64,1})
    p,u,v = SpeedPressureFromCons(obj,conservatives);
    flux = zeros(4,2);
    flux[1,1]=conservatives[2];
    flux[2,1]=conservatives[2]*u+p;
    flux[3,1]=conservatives[2]*v;
    flux[4,1]=(conservatives[4]+p)*u;
    flux[1,2]=conservatives[3];
    flux[2,2]=conservatives[3]*u;
    flux[3,2]=conservatives[3]*v+p;
    flux[4,2]=(conservatives[4]+p)*v;
    return flux;
end

function getCFL(obj::Euler,state::Array{Float64,1},dt,dx)
    rho,p,u,v,a,T = primFromCons(obj,state);
    speedMinU = u-a;
    speedMaxU = u+a;
    speedMinV = v-a;
    speedMaxV = v+a;
    speedU=max(abs(speedMinU),abs(speedMaxU));
    speedV=max(abs(speedMinV),abs(speedMaxV));
    speed = max(abs(speedU),abs(speedV))
    cfl=speed*dt/dx; 
    return cfl;
end