clear all;

Goomrange = 0:6;
muoomrange = -5:3;
G = 433%10.^(linspace(min(Goomrange),max(Goomrange),100*numel(Goomrange)))';
mu = 40*0.894*1e-3;%10.^(linspace(min(muoomrange),max(muoomrange),100*numel(muoomrange)));

alpha = 0.0606;%0.1009;
Rmax = 4.185128557872455e-04;%438.25E-6;

rho = 998.2;
p_inf = 101325;
pR = 92;%10E6/p_inf;
k = 1.4;
c = 1566;
gamma = 0.056;
c0 = sqrt(p_inf/rho);

We = Rmax*p_inf/(2*gamma);
Re = Rmax*rho*c0./mu;
Ca = p_inf./G;%rho*c0^2./G;
C = c/c0;

% nT1 = (1+1.0717E-6/Rmax)*(4.2/(alpha*Rmax)^4.2);
% nT2 = -1.0717E-6/(alpha*Rmax^2);
% nT3 = 3.8274E-5/(alpha*Rmax)^3;
% dT1 = 1;
% dT2 = 6.1271E-3/(alpha*Rmax);
nT1 = (1+1./We)*(3*k)/(alpha)^(3*k);
nT2 = -1./(We*alpha);
nT3 = 4./(Ca*(alpha)^3);
dT1 = 1;
dT2 = 4./(Re*C*alpha);

A = 1/(alpha);

omega_n_star = A*sqrt((nT1+nT2+nT3)*(1./(dT1+dT2)));
omega_n = omega_n_star*c/Rmax;

% beta_t1 = 3.1743E3*alpha*Rmax*omega_n.^2;
% beta_t2 = 2./(alpha*Rmax)./(4.145*alpha*Rmax./mu+2.5397E2);
% beta_tot = sqrt(beta_t1+repmat(beta_t2,size(G)));
beta_t1 = (alpha)*omega_n_star.^2/(2*C);
beta_t2 = 2/(alpha)./(Re*alpha+4/C);
% beta_tot = beta_t1+beta_t2;
beta_tot_star = (beta_t1+repmat(beta_t2,size(G)));

omega_star = sqrt(omega_n_star.^2-beta_tot_star.^2);

%%
figure(1)
imagesc(log10(omega_n_star)); axis image; colorbar; set(gca,'YDir','normal'); ylabel('log(G [Pa])'); xlabel('log(\mu [Pa s])'); 
title('log(\omega_N* [])');
set(gca,'Ytick',linspace(1,length(G),numel(Goomrange))); set(gca,'YTickLabel',Goomrange);
set(gca,'Xtick',linspace(1,length(mu),numel(muoomrange))); set(gca,'XTickLabel',muoomrange);

figure(2)
imagesc(log10(beta_tot_star)); axis image; colorbar; set(gca,'YDir','normal');ylabel('log(G [Pa])'); xlabel('log(\mu [Pa s])');
title('log(\beta_{tot}* [])');
set(gca,'Ytick',linspace(1,length(G),numel(Goomrange))); set(gca,'YTickLabel',Goomrange);
set(gca,'Xtick',linspace(1,length(mu),numel(muoomrange))); set(gca,'XTickLabel',muoomrange);

figure(3)
imagesc(log10(real(omega_star))); axis image; colorbar; set(gca,'YDir','normal'); ylabel('log(G [Pa])'); xlabel('log(\mu [Pa s])');
title('log(\omega* [])');
set(gca,'Ytick',linspace(1,length(G),numel(Goomrange))); set(gca,'YTickLabel',Goomrange);
set(gca,'Xtick',linspace(1,length(mu),numel(muoomrange))); set(gca,'XTickLabel',muoomrange);

% 
% t = linspace(0,100*3.7E-6,1000);
% tstar = t*c0/Rmax;
% idx = [100,400];
% RR = alpha+(1-alpha)*(cos(omega)*t)*exp(-beta_tot*t);
% plot(tstar,RR);
% 
%tc = real(1./(beta_tot_star-(1-alpha)*omega_n_star.^2*C/pR)+pi./omega_star);
% 
% figure(4)
% tc(tc<0) = NaN;
% imagesc(log10(tc)); axis image; colorbar; set(gca,'YDir','normal'); ylabel('log(G [Pa])'); xlabel('log(\mu [Pa s])');
% title('log(tc [])');
% set(gca,'Ytick',linspace(1,length(G),numel(Goomrange))); set(gca,'YTickLabel',Goomrange);
% set(gca,'Xtick',linspace(1,length(mu),numel(muoomrange))); set(gca,'XTickLabel',muoomrange);
% 

%% 
tstar=linspace(0,1.4E-3,1001);
Rstar = 0.5*cos(omega_n_star*tstar).*exp(-beta_tot_star*tstar)+0.5;
plot(tstar,Rstar); set(gca,'YLim',[0 1]);
hold on;

