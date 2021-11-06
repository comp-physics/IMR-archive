figure; plot(t,Rnew(expts,:)*1e-6,'rs')
hold on; plot(t0+t2_num,R2_num)

title(['G=',num2str(G),'; \mu=',num2str(mu),'; \lambda_\nu=',num2str(lambda_nu), ...
    '; \alpha=',num2str(alpha)],'FontWeight','Normal');
set(gca,'fontsize',18)

% title(['G=',num2str(10^G_ooms(i)),' \mu=',num2str(0.014),', \alpha=',num2str(0)],'FontWeight','Normal');
% set(gca,'fontsize',18)

%%
% fp = './data/1pt3kPa_PA/'; cd([fp num2str(expt)]); 
% load(['Rfit_stiff_exp',num2str(expt),'.mat'],'Rfit_pp','t2Start');
fp = ['/Users/yangjin/Documents/MATLAB/Franck/IMR-master/data/numsim/',num2str(simNo),'/'];
filename_temp = ['Rfit_numsim',num2str(simNo),'_100.mat'];
load([fp filename_temp]);
% load(['Rfit_soft_exp',num2str(expt),'.mat'],'Rfit_pp','t2Start'); 
% load(['Rfit_collagen1_exp',num2str(expt),'.mat'],'Rfit_pp','t2Start'); 
tNew =  0:1e-7:t(end); Rfit_exp = ppval(Rfit_pp,tNew); % figure; plot(tNew,Rfit,'rs');

% ====== comment codes ======
% h_exp = envelope(Rfit,30,'peaks');
% hold on; plot(tNew,h_exp);
% h_exp = envelope(t,Rnew(expts,:)*1e-6,'linear');
% hold on; plot(t,h_exp);
d_exp =  Rfit_exp; 
Ts = tNew(end)/length(d_exp); Fs = 1/Ts; Fn = Fs/2; tv_exp = (0:length(d_exp)-1)*Ts; % figure; plot(tv,d,'rs');
st = find(d_exp > 1e-5, 1, 'first'); tvst_exp = tv_exp(st);
d_exp = d_exp(st:end); tv_exp = tv_exp(st:end); L = length(d_exp);
[pks, locs] = findpeaks( (d_exp), tv_exp, 'MinPeakDist',1E-6);
peakVal_exp = pks(1:10)';
peakLoc_exp = locs(1:10)'+tvst_exp ;
q = [locs' ones(size(locs'))]\log(abs(pks))';          % Initial Parameter Estimaes
fitfcn = @(b,t) b(1) .* exp(b(2).*t) + b(3);           % Fit Function
SSECF = @(b) sum((pks' - fitfcn(b,locs')).^2);
[B_exp,SSE] = fminsearch(SSECF, [pks(1); q(1); 0]);
figure; % plot(t,Rnew(expts,:)*1e-6,'rs'); hold on; 
plot(tv_exp+tvst_exp , fitfcn(B_exp,tv_exp),'r--'); hold on;
plot(peakLoc_exp,peakVal_exp,'rd') ;hold on;
plot(tv_exp+tvst_exp,d_exp,'r'); hold off; 
grid off; grid; set(gca,'fontsize',18); grid minor; 
% title(['G=',num2str(10^G_ooms),' \mu=',num2str(10^mu_ooms),', \alpha=',num2str(alpha)],'FontWeight','Normal');
title('Experiment data','FontWeight','Normal');
set(gca,'fontsize',18)

[pksDown, locsDown] = findpeaks( (-d_exp), tv_exp, 'MinPeakDist',1E-6);
peakDownVal_exp = pksDown(1:10-1)';
peakDownLoc_exp = locsDown(1:10-1)'+tvst_exp ;

%%
deltat = 1e-7; t3 = t0:deltat:tspan+t0;
Rfit_pp_num = interp1(t0+t2_num,R2_num,'pchip','pp'); Rfit_num = ppval(Rfit_pp_num,t3);
d = Rfit_num;
Ts = (t3(end)-t0)/length(d); Fs = 1/Ts; Fn = Fs/2; tv_num = (0:length(d)-1)*Ts; 
st = find(d > 1e-5, 1, 'first'); tvst_num = tv_num(st);
d = d(st:end); tv_num = tv_num(st:end); L = length(d); % figure; plot(tv_num,d);
[pks, locs] = findpeaks( (d), tv_num, 'MinPeakDist',1e-6);
pks = [d(1) pks]; locs = [0 locs]; % hold on; plot(locs,pks);
peakVal_num = pks(1:min([5,length(pks)]))';
peakLoc_num = locs(1:min([5,length(pks)]))'+tvst_num ;
% q = [locs' ones(size(locs'))]\log(abs(pks))';           % Initial Parameter Estimaes
% fitfcn = @(b,t) b(1) .* exp(b(2).*t) + b(3);            % Fit Function
% SSECF = @(b) sum((pks' - fitfcn(b,locs')).^2);          % Cost Function
% [B_num,SSE] = fminsearch(SSECF, [pks(1); q(1); 0]); 

figure; plot(t0+tv_num+tvst_num, d,'b'); hold on; % plot(t0+tv+tvst, fitfcn(B_num,tv)); hold off; grid
plot(tv_exp+tvst_exp,d_exp,'r'); hold on;
plot(t0+tvst_num+peakLoc_num,peakVal_num,'bo') ;hold on;
plot(tv_exp+tvst_num , fitfcn(B_exp,tv_exp),'r--'); hold on;
plot(peakLoc_exp,peakVal_exp,'rd') ;hold on;
grid off; grid; set(gca,'fontsize',18); grid minor; 
% title(['G=',num2str(10^G_ooms(i)),' \mu=',num2str(10^mu_ooms(j)),', \alpha=',num2str(alpha)],'FontWeight','Normal');
% title('Experiment data','FontWeight','Normal');
title(['G=',num2str(G),'; \mu=',num2str(mu),'; \lambda_\nu=',num2str(lambda_nu), ...
    '; \alpha=',num2str(alpha)],'FontWeight','Normal');
set(gca,'fontsize',18); 
% axis([0,4e-4,0,2.5e-4]); 
axis([0,2.5e-4,0,2.5e-4]);
lgd = legend('Fung model Kelvin-Voigt','Experiment data');
% 'Neo-Hookean Kelvin-Voigt'

%%
% [peakVal_exp peakLoc_exp]
% [peakVal_num peakLoc_num+(t0+tvst_num)]

%% Make movie

% workingDir = '/Users/yangjin/Documents/MATLAB/Franck/IMR-master';
% % outputVideo = VideoWriter(fullfile(workingDir,'collagen1_expt6_Fung_ln_change_pkloc'));
% % outputVideo = VideoWriter(fullfile(workingDir,'stiff_expt3_Fung_ln_change_pkloc'));
% 
% outputVideo.FrameRate = 1;
% open(outputVideo);
% 
% imgindex = [1:1:5];
% for tempi = imgindex
% figure(tempi); fig=gcf; frame=getframe(gcf); writeVideo(outputVideo,frame);
% end
% 
% close(outputVideo);

