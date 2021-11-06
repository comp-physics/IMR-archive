
% if varagin(t2_num)  
% else
    t2=t2_num; % R2=R2_num;
% end

%% Plot C mass fraction of gas 
figure; plot(t2,C(:,1)); hold on; plot(t2,C(:,400)); hold on; plot(t2,C(:,450));
hold on; plot(t2,C(:,475)); hold on; plot(t2,C(:,500)); set(gca,'yscale','log'); set(gca,'fontsize',24)
% title('IMRsolver1: C-t','FontWeight','normal')
title('IMRsolver2: C-t','FontWeight','normal')
% title('Exact: C-t','FontWeight','normal')
% axis([0,tspan,1e-3,1e1]); 
lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot T temperature inside bubble
figure; plot(t2,T(:,1)); hold on; plot(t2,T(:,400));hold on; plot(t2,T(:,450))
hold on; plot(t2,T(:,475));hold on; plot(t2,T(:,500));set(gca,'yscale','log');
set(gca,'fontsize',24);
% title('IMRsolver1: T-t','FontWeight','normal')
title('IMRsolver2: T-t','FontWeight','normal')
% title('Exact: T-t','FontWeight','normal')
lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot P pressure inside bubble
 figure; t2=t2_num; plot(t2,P)
set(gca,'yscale','log')
% title('IMRsolver1: P-t','FontWeight','normal')
title('IMRsolver2: P-t','FontWeight','normal')
% title('Exact: P-t','FontWeight','normal')
set(gca,'fontsize',24); axis([0,2.5e-4,0,1e8])

% lgd = legend('IMR Solver New','IMR Solver Old','IMR Exact')

%% Plot R
figure; plot(t2,R2)
% title('IMRsolver1: R-t','FontWeight','normal')
title('IMRsolver2: R-t','FontWeight','normal')
% title('Exact: R-t','FontWeight','normal')
set(gca,'fontsize',24); 
% axis([0,2.5e-4,0,2.5e-4])

%% Plot Tm medium temperature outside bubble
% Nm = NTM-1;
% deltaYm = -2/Nm;
% j = 1:1:Nm+1;
% xk = (1+(j-1)*deltaYm)';
% L=2; yk2 = ((2./(xk+1)-1)*L+1);
% 
% figure; plot(t2,Tm(:,1)); hold on; plot(t2,Tm(:,2)); hold on; plot(t2,Tm(:,3)); hold on;
% plot(t2,Tm(:,4)); hold on; plot(t2,Tm(:,5)); hold on; plot(t2,Tm(:,6)); hold on; 
% plot(t2,Tm(:,7)); hold on; plot(t2,Tm(:,8)); hold on; plot(t2,Tm(:,9)); hold on; plot(t2,Tm(:,10));
% set(gca,'yscale','log')
% set(gca,'fontsize',24)
% title('Jon-IMRsolver: Tm/T{\infty}-t','FontWeight','normal')
% lgd=legend('All yk');



%%  ************** IMR solver Old *************
%% Plot C mass fraction of gas 
  t2=t2_num; R2=R2_num;figure; plot(t2,C(:,1)); hold on; plot(t2,C(:,400)); hold on; plot(t2,C(:,450));
hold on; plot(t2,C(:,475)); hold on; plot(t2,C(:,500)); set(gca,'yscale','log'); set(gca,'fontsize',24)
  title('IMRsolver1: C-t','FontWeight','normal')
axis([0,tspan,1e-3,1e1]); lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot T temperature inside bubble
  t2=t2_num; R2=R2_num; figure; plot(t2,T(:,1)); hold on; plot(t2,T(:,400));hold on; plot(t2,T(:,450))
hold on; plot(t2,T(:,475));hold on; plot(t2,T(:,500));set(gca,'yscale','log');set(gca,'fontsize',24);
title('IMRsolver1: T-t','FontWeight','normal')
lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot P pressure inside bubble
  figure; t2=t2_num; R2=R2_num; hold on; plot(t2,P)
set(gca,'yscale','log')
  title('IMRsolver1: P-t','FontWeight','normal')
set(gca,'fontsize',24); axis([0,2.5e-4,0,1e8])

% lgd = legend('IMR Solver New','IMR Solver Old','IMR Exact')
% title('G=0; \mu=0.001; \alpha=0; \lambda_\nu=1')
% fig_NHKV_G2970_mu5_a5_ln1_P-t_comp

%% Plot R
figure; plot(t2,R2)
  title('IMRsolver1: R-t','FontWeight','normal')
set(gca,'fontsize',24); axis([0,2.5e-4,0,2.5e-4])


 


%%  ************** IMR Exact *************
%% Plot C mass fraction of gas 
  t2=t2_num; R2=R2_num;figure; plot(t2,C(:,1)); hold on; plot(t2,C(:,400)); hold on; plot(t2,C(:,450));
hold on; plot(t2,C(:,475)); hold on; plot(t2,C(:,500)); set(gca,'yscale','log'); set(gca,'fontsize',24)
  title('Exact: C-t','FontWeight','normal')
axis([0,tspan,1e-3,1e1]); lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot T temperature inside bubble
  t2=t2_num; R2=R2_num; figure; plot(t2,T(:,1)); hold on; plot(t2,T(:,400));hold on; plot(t2,T(:,450))
hold on; plot(t2,T(:,475));hold on; plot(t2,T(:,500));set(gca,'yscale','log');set(gca,'fontsize',24);
title('Exact: T-t','FontWeight','normal')
lgd=legend('R*=0','R*=0.8','R*=0.9','R*=0.95','R*=1');

%% Plot P pressure inside bubble
 figure;  hold on;   plot(t2,P)
set(gca,'yscale','log')
  title('Exact: P-t','FontWeight','normal')
set(gca,'fontsize',24); axis([0,2.5e-4,0,1e8])

% lgd = legend('IMR Solver New','IMR Solver Old','IMR Exact');

%% Plot R
figure; plot(t2,R2)
  title('Exact: R-t','FontWeight','normal')
set(gca,'fontsize',24); axis([0,2.5e-4,0,2.5e-4])





