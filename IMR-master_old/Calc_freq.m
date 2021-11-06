function omega = Calc_freq(G,mu,rho,R0,k,P_inf,T_inf,gamma)

Pv = Pvsat(1*T_inf); 
term1 = 1/(rho*R0^2)*(3*k*(P_inf-Pv)-(3*k-1)*2*gamma/(rho*R0));
term2 = -8*mu^2/(rho^2*R0^4);
term3 = 4*G/(rho*R0^2);

omega=sqrt(term1+term2+term3);

end


