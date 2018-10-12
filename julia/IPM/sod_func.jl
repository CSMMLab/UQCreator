function sod_func(P,rho_l,P_l,u_l,rho_r,P_r,u_r,gamma)
#defines function to be used in analytic_sod
#Initial conditions


#z = (P / P_r - 1.0);
#c1 = sqrt(gamma * P_l / rho_l);
#c5 = sqrt(gamma * P_r / rho_r);

#gm1 = gamma - 1.0;
#gp1 = gamma + 1.0;
#g2 = 2.0 * gamma;

#fact = gm1 / g2 * (c5 / c1) * z / sqrt(1.0 + gp1 / g2 * z)
#fact = (1.0 - fact)^(g2 / gm1)

#return P_l * fact - P_r


mu = sqrt( (gamma-1.0)/(gamma+1.0) );

#y = (P - P_r)*(( ((1 - mu^2)^2)*((rho_r*(P + mu*mu*P_r))^-1) )^(0.5)) - 2*(sqrt(gamma)/(gamma - 1))*(1 - (P)^( (gamma - 1)/(2*gamma)));
y = (P - P_r)*((1.0-mu*mu)*(rho_r*(P + mu*mu*P_r))^(-1))^(0.5)-(P_l^( (gamma-1.0)/(2*gamma))-P^( (gamma-1.0)/(2*gamma)))*(((1.0-mu*mu*mu*mu)*P_l^(1.0/gamma))*(mu*mu*mu*mu*rho_l)^(-1))^(0.5);


return y

end
