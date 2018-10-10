using Roots
include("sod_func.jl")
function analytic_sod(t,x0,rho_l,P_l,u_l,rho_r,P_r,u_r,gamma,x)
#to solve Sod's Shock Tube problem
#reference: "http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html"
#   |       |   |     |         |
#   |       |   |     |         |
#   |       |   |     |         |
#___|_______|___|_____|_________|_______________
#   x1      x2  x0    x3        x4
#
#input require: t (time)

mu = sqrt( (gamma-1)/(gamma+1) );

#speed of sound
c_l = ( (gamma*P_l/rho_l))^(0.5);
c_r = ( (gamma*P_r/rho_r))^(0.5);

#P_post = find_zero(xi->sod_func(xi,rho_l,P_l,u_l,rho_r,P_r,u_r,gamma), (0.0,10^7), Bisection())
P_post = find_zero(xi->sod_func(xi,rho_l,P_l,u_l,rho_r,P_r,u_r,gamma), (P_l,P_r), Bisection())
rho_post = rho_r*(( (P_post/P_r) + mu^2 )/(1 + mu*mu*(P_post/P_r)));
rho_middle = (rho_l)*((P_post/P_l))^(1/gamma);



gm1 = gamma - 1.0;
gp1 = gamma + 1.0;
gmfac1 = 0.5 * gm1 / gamma;
gmfac2 = 0.5 * gp1 / gamma;

z = (P_post / P_r - 1.0)
c5 = sqrt(gamma * P_r / rho_r)

fact = sqrt(1.0 + gmfac2 * z);

# shock speed
w = c5 * fact;

u4 = c5 * z / (gamma * fact);

#Key Positions
x4 = x0 + w * t



x1 = x0 - c_l*t;
x3 = x0 + u4*t;
#determining x2
c_2 = c_l - ((gamma - 1.0)/2)*u4;
x2 = x0 + (u4 - c_2)*t;

#start setting values
n_points = length(x);

data_x = x';
data_rho = zeros(n_points);   #density
data_P = zeros(n_points); #pressure
data_u = zeros(n_points); #velocity
data_e = zeros(n_points); #internal energy

for index = 1:n_points
    if data_x[index] < x1
        #Solution b4 x1
        data_rho[index] = rho_l;
        data_P[index] = P_l;
        data_u[index] = u_l;
    elseif (x1 <= data_x[index] && data_x[index] <= x2)
        #Solution b/w x1 and x2
        c = mu*mu*((x0 - data_x[index])/t) + (1 - mu*mu)*c_l; 
        data_rho[index] = rho_l*(c/c_l)^(2/(gamma - 1));
        data_P[index] = P_l*((data_rho[index]/rho_l))^gamma;
        data_u[index] = (1 - mu*mu)*( (-(x0-data_x[index])/t) + c_l);
    elseif (x2 <= data_x[index] && data_x[index] <= x3)
        #Solution b/w x2 and x3
        data_rho[index] = rho_middle;
        data_P[index] = P_post;
        data_u[index] = u4;
    elseif (x3 <= data_x[index] && data_x[index] <= x4)
        #Solution b/w x3 and x4
        data_rho[index] = rho_post;
        data_P[index] = P_post;
        data_u[index] = u4;
    elseif x4 < data_x[index]
        #Solution after x4
        data_rho[index] = rho_r;
        data_P[index] = P_r;
        data_u[index] = u_r;
    end
    data_e[index] = data_P[index]/((gamma - 1)*data_rho[index]);
end
return data_rho,data_u,data_P,data_e;
end
